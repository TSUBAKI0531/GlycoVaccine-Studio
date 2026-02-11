import streamlit as st
import pandas as pd
import json
from io import BytesIO
import streamlit.components.v1 as components
from your_module import (
    GlycoConjugateWorkflow, 
    AntibodyDesigner, 
    AntibodyDockingWorkflow, 
    LightweightHotSpotAnalyzer
)

# --- Session State 初期化 ---
if 'analysis_history' not in st.session_state: st.session_state.analysis_history = []
if 'last_antigen_prot' not in st.session_state: st.session_state.last_antigen_prot = ""

# --- 3D可視化用関数 ---
def show_3d_model(cif_text):
    safe_cif = cif_text.replace("`", "\\`").replace("$", "\\$")
    html_code = f"""
    <div id="container" style="height: 500px; width: 100%;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        $(function() {{
            let viewer = $3Dmol.createViewer($("#container"), {{backgroundColor: "white"}});
            viewer.addModel(`{safe_cif}`, "mcif");
            viewer.setStyle({{cartoon: {{color: 'spectrum'}}}});
            viewer.setStyle({{hetflag: true}}, {{stick: {{radius: 0.3}}}});
            viewer.zoomTo(); viewer.render();
        }});
    </script>"""
    components.html(html_code, height=520)

# --- アプリ設定 ---
st.set_page_config(page_title="GlycoVaccine Studio", layout="wide", page_icon="🧪")
st.title("🧪 GlycoVaccine Studio v1.5")

with st.sidebar:
    st.header("⚙️ Settings")
    job_name = st.text_input("Job Name", "TACA_Project")
    h_id = st.text_input("H-Chain ID", "H")
    l_id = st.text_input("L-Chain ID", "L")

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "🧬 抗原デザイン", "📊 モデル解析", "🛡️ 抗体解析", "🔥 Hot Spot", "🎨 抗体デザイン", "📋 履歴"
])

# --- Tab 1: 抗原デザイン ---
with tab1:
    prot_seq = st.text_area("Antigen Protein Sequence")
    smiles = st.text_input("Glycan SMILES")
    bond_idx = st.number_input("Bonding Residue Index", value=1)
    if st.button("Save Antigen Info"):
        st.session_state.last_antigen_prot = prot_seq
        st.session_state.last_smiles = smiles
        st.session_state.last_bond_idx = bond_idx
        st.success("抗原情報を保存しました。")

# --- Tab 2: モデル解析 ---
with tab2:
    st.header("📊 Model Evaluation (SASA)")
    uploaded = st.file_uploader("Upload CIF models", accept_multiple_files=True, key="eval")
    if uploaded:
        wf = GlycoConjugateWorkflow(job_name)
        for f in uploaded:
            content = f.read().decode("utf-8")
            with open(f.name, "w") as tmp: tmp.write(content)
            res = wf._calculate_sasa(f.name)
            st.write(f"Model: {f.name}, SASA: {res['glycan_sasa']:.2f}")
            st.session_state.analysis_history.append({"Job": job_name, "Type": "SASA", "Target": f.name, "Score": res['glycan_sasa']})

# --- Tab 3: 抗体解析 (復旧箇所) ---
with tab3:
    st.header("🛡️ Paratope Analysis")
    complex_file = st.file_uploader("Upload Complex CIF", key="comp")
    if complex_file:
        content = complex_file.read().decode("utf-8")
        with open("temp_comp.cif", "w") as f: f.write(content)
        
        if st.button("Analyze Paratope"):
            adw = AntibodyDockingWorkflow(job_name, h_chain=h_id, l_chain=l_id)
            df_para = adw.analyze_paratope("temp_comp.cif")
            st.dataframe(df_para)
            show_3d_model(content)

# --- Tab 4: Hot Spot解析 (復旧箇所) ---
with tab4:
    st.header("🔥 Hot Spot Prediction")
    hs_file = st.file_uploader("Upload Structure", key="hs")
    if hs_file:
        content = hs_file.read().decode("utf-8")
        with open("temp_hs.cif", "w") as f: f.write(content)
        if st.button("Run Hot Spot Scan"):
            analyzer = LightweightHotSpotAnalyzer("temp_hs.cif", h_chain=h_id, l_chain=l_id)
            res = analyzer.run_contact_density_scan()
            st.dataframe(res.head(10))
            st.bar_chart(res.set_index("Residue")["HotSpot_Score"])

# --- Tab 5: 抗体デザイン ---
with tab5:
    st.header("🎨 Antibody Binder Design")
    designer = AntibodyDesigner()
    motif = st.selectbox("Target Motif", ["Tn Antigen"])
    ranked = designer.get_ranked_candidates(motif)
    st.table(pd.DataFrame(ranked)[["name", "Score"]])
    
    selected = st.selectbox("Select Candidate", [r["name"] for r in ranked])
    best = next(r for r in ranked if r["name"] == selected)
    
    if st.button("Generate AF3 Full JSON"):
        if st.session_state.last_antigen_prot:
            wf = GlycoConjugateWorkflow(job_name)
            full_json = wf.create_full_complex_json(
                job_name, st.session_state.last_antigen_prot, 
                st.session_state.last_smiles, st.session_state.last_bond_idx, 
                best["H_AA"], best["L_AA"]
            )
            st.download_button("Download Full JSON", json.dumps(full_json, indent=2), f"{job_name}_full.json")
        else:
            st.warning("先にタブ1で抗原情報を保存してください。")

# --- Tab 6: 履歴 ---
with tab6:
    if st.session_state.analysis_history:
        df_h = pd.DataFrame(st.session_state.analysis_history)
        st.dataframe(df_h)
        output = BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            df_h.to_excel(writer, index=False)
        st.download_button("Download Excel", output.getvalue(), "history.xlsx")