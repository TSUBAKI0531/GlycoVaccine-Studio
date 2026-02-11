import streamlit as st
import pandas as pd
import json
from your_module import GlycoConjugateWorkflow, AntibodyDesigner

# --- Session State ---
if 'last_antigen_prot' not in st.session_state: st.session_state.last_antigen_prot = ""

st.set_page_config(page_title="GlycoVaccine Studio", layout="wide", page_icon="🧪")
st.title("🧪 GlycoVaccine Studio v1.5")

with st.sidebar:
    st.header("⚙️ Settings")
    platform = st.radio("Target Platform", ["AlphaFold Server (Web)", "AlphaFold 3 (Standalone)"])
    job_name = st.text_input("Job Name", "TACA_Project")

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["🧬 抗原デザイン", "📊 モデル解析", "🛡️ 抗体解析", "🔥 Hot Spot", "🎨 抗体デザイン", "📋 履歴"])

with tab1:
    st.header("🧬 Antigen Design")
    prot_seq = st.text_area("Antigen Protein Sequence (Carrier)")
    
    if platform == "AlphaFold Server (Web)":
        glycan_input = st.text_input("Glycan CCD Code", help="例: A2G (GalNAc), NAG (GlcNAc)")
    else:
        glycan_input = st.text_input("Glycan SMILES")
    
    if st.button("Save Antigen Info"):
        st.session_state.last_antigen_prot = prot_seq
        st.session_state.last_glycan = glycan_input
        st.success("抗原情報を保存しました。")

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
            st.download_button("Download JSON", json.dumps(full_json, indent=2), f"{job_name}_{mode}.json")
        else:
            st.warning("Tab 1で情報を保存してください。")