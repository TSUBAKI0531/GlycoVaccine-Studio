import streamlit as st
import pandas as pd
import io
import os
from io import BytesIO
import streamlit.components.v1 as components
from your_module import (
    GlycoConjugateWorkflow, 
    AntibodyDockingWorkflow, 
    LightweightHotSpotAnalyzer
)

# --- 履歴保持の初期化 (Session State) ---
if 'analysis_history' not in st.session_state:
    st.session_state.analysis_history = []

# --- 3D可視化用の補助関数 ---
def show_3d_model(cif_text):
    safe_cif = cif_text.replace("`", "\\`").replace("$", "\\$")
    html_code = f"""
    <div id="container" style="height: 500px; width: 100%; position: relative;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        $(function() {{
            let viewer = $3Dmol.createViewer($("#container"), {{backgroundColor: "white"}});
            let data = `{safe_cif}`;
            viewer.addModel(data, "mcif");
            viewer.setStyle({{cartoon: {{color: 'spectrum'}}}});
            viewer.setStyle({{hetflag: true}}, {{stick: {{colorscheme: 'magentaCarbon', radius: 0.3}}}});
            viewer.zoomTo();
            viewer.render();
        }});
    </script>
    """
    components.html(html_code, height=520)

# --- Streamlit設定 ---
st.set_page_config(page_title="GlycoVaccine Studio", page_icon="🧪", layout="wide")
st.title("🧪 GlycoVaccine Studio")

with st.sidebar:
    st.header("⚙️ Project Settings")
    job_name = st.text_input("Job Name", "TACA_Project_V1")
    h_chain_global = st.text_input("Antibody Heavy Chain ID", "H")
    l_chain_global = st.text_input("Antibody Light Chain ID", "L")

# タブに「履歴＆出力」を追加
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🧬 抗原デザイン", "📊 モデル解析", "🛡️ 抗体ドッキング解析", "🔥 Hot Spot解析", "📋 履歴＆出力"
])

# --- Tab 1: 抗原デザイン (変更なし) ---
with tab1:
    st.header("🧬 Antigen Design")
    col1, col2 = st.columns(2)
    with col1:
        prot_seq = st.text_area("Carrier Protein Sequence")
        smiles = st.text_input("Glycan-Linker SMILES")
    with col2:
        smarts = st.text_input("Terminal SMARTS", "C(=O)N")
        bond_idx = st.number_input("Bonding Residue Index", value=1)
    
    if st.button("Generate AF3 JSON"):
        if prot_seq and smiles:
            wf = GlycoConjugateWorkflow(job_name)
            json_path = wf.prepare_af3_input(prot_seq, smiles, bond_idx, smarts)
            with open(json_path, "r") as f:
                st.download_button("Download JSON", f.read(), f"{job_name}.json")
            st.success("JSON Created!")

# --- Tab 2: モデル解析 ---
with tab2:
    st.header("📊 Model Evaluation")
    uploaded_models = st.file_uploader("Upload CIF models", accept_multiple_files=True, key="eval")
    
    if uploaded_models:
        results = []
        wf = GlycoConjugateWorkflow(job_name)
        for file in uploaded_models:
            content = file.read().decode("utf-8")
            with open(file.name, "w") as f: f.write(content)
            sasa = wf._calculate_sasa(file.name)
            contacts = wf.analyze_interactions(file.name)
            score = sasa["glycan_sasa"] / (len(contacts) + 1)
            results.append({"Model": file.name, "Exposure_Score": score, "Content": content})
        
        df_eval = pd.DataFrame(results).sort_values("Exposure_Score", ascending=False)
        st.dataframe(df_eval[["Model", "Exposure_Score"]])
        
        if st.button("結果を履歴に追加", key="save_eval"):
            for _, row in df_eval.iterrows():
                st.session_state.analysis_history.append({
                    "Job": job_name, "Type": "Model Evaluation", 
                    "Target": row["Model"], "Score": row["Exposure_Score"]
                })
            st.toast("履歴に保存しました！")

        if st.checkbox("Show Best Model 3D View"):
            show_3d_model(df_eval.iloc[0]["Content"])

# --- Tab 3: 抗体解析 ---
with tab3:
    st.header("🛡️ Paratope Analysis")
    complex_file = st.file_uploader("Upload Complex CIF", key="comp")
    if complex_file:
        content = complex_file.read().decode("utf-8")
        with open("temp.cif", "w") as f: f.write(content)
        
        if st.button("Analyze Paratope"):
            adw = AntibodyDockingWorkflow(job_name, h_chain=h_chain_global, l_chain=l_chain_global)
            df_para = adw.analyze_paratope("temp.cif")
            st.dataframe(df_para)
            
            # 履歴に追加
            st.session_state.analysis_history.append({
                "Job": job_name, "Type": "Paratope Analysis", 
                "Target": complex_file.name, "Score": len(df_para) # 接触残基数をスコアとして記録
            })
            show_3d_model(content)

# --- Tab 4: Hot Spot解析 ---
with tab4:
    st.header("🔥 Hot Spot Prediction")
    hs_file = st.file_uploader("Upload Structure", key="hs")
    if hs_file:
        if st.button("Run Scan"):
            analyzer = LightweightHotSpotAnalyzer(hs_file, h_chain=h_chain_global, l_chain=l_chain_global)
            res = analyzer.run_contact_density_scan()
            st.dataframe(res.head(10))
            st.bar_chart(res.set_index("Residue")["HotSpot_Score"])
            
            # 履歴に追加
            top_res = res.iloc[0]
            st.session_state.analysis_history.append({
                "Job": job_name, "Type": "Hot Spot", 
                "Target": hs_file.name, "Score": top_res["HotSpot_Score"]
            })

# --- Tab 5: 履歴＆出力 (新設) ---
with tab5:
    st.header("📋 Analysis History & Export")
    
    if st.session_state.analysis_history:
        df_history = pd.DataFrame(st.session_state.analysis_history)
        st.dataframe(df_history, use_container_width=True)
        
        col1, col2, col3 = st.columns(3)
        
        # CSV出力
        with col1:
            csv = df_history.to_csv(index=False).encode('utf-8')
            st.download_button("Download CSV", csv, "history.csv", "text/csv")
            
        # Excel出力
        with col2:
            output = BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df_history.to_excel(writer, index=False, sheet_name='Summary')
            st.download_button("Download Excel", output.getvalue(), "history.xlsx")
            
        # 履歴削除
        with col3:
            if st.button("履歴をリセット", color="primary"):
                st.session_state.analysis_history = []
                st.rerun()
    else:
        st.info("解析履歴はまだありません。各タブで解析を実行してください。")