import streamlit as st
import pandas as pd
import json
import streamlit.components.v1 as components
from your_module import (
    GlycoConjugateWorkflow, AntibodyDesigner, 
    AntibodyDockingWorkflow, LightweightHotSpotAnalyzer
)

# --- Session State ---
if 'analysis_history' not in st.session_state: st.session_state.analysis_history = []
if 'last_antigen_prot' not in st.session_state: st.session_state.last_antigen_prot = ""
if 'last_glycan' not in st.session_state: st.session_state.last_glycan = ""

# --- 3D 可視化用コンポーネント ---
def show_3d_viewer(cif_text):
    safe_cif = cif_text.replace("`", "\\`").replace("$", "\\$")
    html_code = f"""
    <div id="container" style="height: 400px; width: 100%;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        $(function() {{
            let viewer = $3Dmol.createViewer($("#container"), {{backgroundColor: "white"}});
            viewer.addModel(`{safe_cif}`, "mcif");
            viewer.setStyle({{cartoon: {{color: 'spectrum'}}}});
            viewer.zoomTo(); viewer.render();
        }});
    </script>"""
    components.html(html_code, height=420)

st.set_page_config(page_title="GlycoVaccine Studio", layout="wide", page_icon="🧪")
st.title("🧪 GlycoVaccine Studio v1.5")

with st.sidebar:
    st.header("⚙️ Settings")
    platform = st.radio("Target Platform", ["AlphaFold Server (Web)", "AlphaFold 3 (Standalone)"])
    job_name = st.text_input("Job Name", "TACA_Project")
    h_id = st.text_input("H-Chain ID", "H")
    l_id = st.text_input("L-Chain ID", "L")

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "🧬 抗原デザイン", "📊 モデル解析", "🛡️ 抗体解析", "🔥 Hot Spot", "🎨 抗体デザイン", "📋 履歴"
])

# --- Tab 1 ---
with tab1:
    st.header("🧬 Antigen Design")
    prot_seq = st.text_area("Antigen Protein Sequence (Carrier)", value=st.session_state.last_antigen_prot)
    glycan_input = st.text_input("Linker + Glycan SMILES", value=st.session_state.last_glycan)
    if st.button("Save Antigen Info"):
        st.session_state.last_antigen_prot = prot_seq
        st.session_state.last_glycan = glycan_input
        st.success("抗原情報を保存しました。")

# --- Tab 2 ---
with tab2:
    st.header("📊 Model Evaluation (SASA)")
    uploaded_eval = st.file_uploader("Upload CIF models for evaluation", accept_multiple_files=True)
    if uploaded_eval:
        wf = GlycoConjugateWorkflow(job_name)
        for f in uploaded_eval:
            content = f.read().decode("utf-8")
            with open("temp.cif", "w") as tmp: tmp.write(content)
            res = wf._calculate_sasa("temp.cif")
            st.metric(label=f"SASA: {f.name}", value=f"{res['glycan_sasa']} Å²")
            st.session_state.analysis_history.append({"Job": job_name, "Target": f.name, "SASA": res['glycan_sasa']})

# --- Tab 3 ---
with tab3:
    st.header("🛡️ Paratope Analysis")
    uploaded_comp = st.file_uploader("Upload Complex CIF", key="comp")
    if uploaded_comp:
        content = uploaded_comp.read().decode("utf-8")
        with open("temp_comp.cif", "w") as f: f.write(content)
        if st.button("Run Paratope Analysis"):
            adw = AntibodyDockingWorkflow(job_name, h_id, l_id)
            df_para = adw.analyze_paratope("temp_comp.cif")
            st.subheader("Potential Interactions")
            st.dataframe(df_para)
            show_3d_viewer(content)

# --- Tab 4 ---
with tab4:
    st.header("🔥 Hot Spot Prediction")
    uploaded_hs = st.file_uploader("Upload Structure for Scan", key="hs")
    if uploaded_hs:
        if st.button("Run Hot Spot Scan"):
            analyzer = LightweightHotSpotAnalyzer(uploaded_hs.name, h_id, l_id)
            df_hs = analyzer.run_contact_density_scan()
            st.subheader("Hot Spot Ranking")
            st.dataframe(df_hs)
            st.bar_chart(df_hs.set_index("Residue")["HotSpot_Score"])

# --- Tab 5 ---
with tab5:
    st.header("🎨 Antibody Design")
    designer = AntibodyDesigner()
    motif = st.selectbox("Target Motif", ["Tn Antigen"])
    ranked = designer.get_ranked_candidates(motif)
    st.table(pd.DataFrame(ranked)[["name", "Score"]])
    selected = st.selectbox("Select Candidate", [r["name"] for r in ranked])
    best = next(r for r in ranked if r["name"] == selected)
    if st.button("Generate AF3 JSON"):
        if st.session_state.last_antigen_prot:
            wf = GlycoConjugateWorkflow(job_name)
            mode = "Web" if platform == "AlphaFold Server (Web)" else "Standalone"
            full_json = wf.create_full_complex_json(
                job_name, st.session_state.last_antigen_prot, 
                st.session_state.last_glycan, best["H_AA"], best["L_AA"], mode=mode
            )
            st.download_button("Download JSON", json.dumps(full_json, indent=2), f"{job_name}_full.json")
        else:
            st.warning("先に「🧬 抗原デザイン」タブで情報を保存してください。")

# --- Tab 6 ---
with tab6:
    if st.session_state.analysis_history:
        st.dataframe(pd.DataFrame(st.session_state.analysis_history))