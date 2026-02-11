import streamlit as st
import pandas as pd
import json
from io import BytesIO
import streamlit.components.v1 as components
from your_module import (
    GlycoConjugateWorkflow, AntibodyDockingWorkflow, 
    LightweightHotSpotAnalyzer, AntibodyDesigner
)

# --- 初期化 ---
if 'analysis_history' not in st.session_state: st.session_state.analysis_history = []
if 'last_antigen_prot' not in st.session_state: st.session_state.last_antigen_prot = ""

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

st.set_page_config(page_title="GlycoVaccine Studio", layout="wide")
st.title("🧪 GlycoVaccine Studio v1.5")

with st.sidebar:
    job_name = st.text_input("Job Name", "TACA_Project")
    h_id = st.text_input("H-Chain ID", "H")
    l_id = st.text_input("L-Chain ID", "L")

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["🧬 抗原デザイン", "📊 モデル解析", "🛡️ 抗体解析", "🔥 Hot Spot", "🎨 抗体デザイン", "📋 履歴"])

with tab1:
    prot_seq = st.text_area("Protein Sequence")
    smiles = st.text_input("Glycan SMILES")
    bond_idx = st.number_input("Bond Index", value=1)
    if st.button("Save & Generate JSON"):
        st.session_state.last_antigen_prot = prot_seq
        st.session_state.last_smiles = smiles
        st.session_state.last_bond_idx = bond_idx
        st.success("抗原情報を保存しました。")

with tab2:
    uploaded = st.file_uploader("Upload CIF", accept_multiple_files=True)
    if uploaded:
        wf = GlycoConjugateWorkflow(job_name)
        for f in uploaded:
            res = wf._calculate_sasa(f.name)
            st.write(f"{f.name}: SASA {res['glycan_sasa']}")
            st.session_state.analysis_history.append({"Job": job_name, "Type": "SASA", "Target": f.name, "Score": res['glycan_sasa']})

with tab5:
    designer = AntibodyDesigner()
    motif = st.selectbox("Target Motif", ["Tn Antigen"])
    ranked = designer.get_ranked_candidates(motif)
    st.table(pd.DataFrame(ranked)[["name", "Score"]])
    
    selected = st.selectbox("Select Candidate", [r["name"] for r in ranked])
    best = next(r for r in ranked if r["name"] == selected)
    
    if st.button("Generate AF3 Full JSON"):
        if st.session_state.last_antigen_prot:
            wf = GlycoConjugateWorkflow(job_name)
            full_json = wf.create_full_complex_json(job_name, st.session_state.last_antigen_prot, st.session_state.last_smiles, st.session_state.last_bond_idx, best["H_AA"], best["L_AA"])
            st.download_button("Download Full JSON", json.dumps(full_json, indent=2), f"{job_name}_full.json")
        else:
            st.warning("先にタブ1で抗原を設計してください。")

with tab6:
    if st.session_state.analysis_history:
        df_h = pd.DataFrame(st.session_state.analysis_history)
        st.dataframe(df_h)
        output = BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            df_h.to_excel(writer, index=False)
        st.download_button("Download Excel", output.getvalue(), "results.xlsx")