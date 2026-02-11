import streamlit as st
import pandas as pd
import json
from io import BytesIO
import streamlit.components.v1 as components
from your_module import GlycoConjugateWorkflow, AntibodyDesigner

if 'analysis_history' not in st.session_state: st.session_state.analysis_history = []
if 'last_antigen_prot' not in st.session_state: st.session_state.last_antigen_prot = ""

st.set_page_config(page_title="GlycoVaccine Studio", layout="wide", page_icon="🧪")
st.title("🧪 GlycoVaccine Studio v1.5")

with st.sidebar:
    st.header("⚙️ Settings")
    job_name = st.text_input("Job Name", "TACA_Project")
    h_id = st.text_input("H-Chain ID", "H")
    l_id = st.text_input("L-Chain ID", "L")

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["🧬 抗原デザイン", "📊 モデル解析", "🛡️ 抗体解析", "🔥 Hot Spot", "🎨 抗体デザイン", "📋 履歴"])

with tab1:
    st.header("🧬 Antigen Design with Linker")
    prot_seq = st.text_area("Antigen Protein Sequence (Carrier)", height=150)
    col1, col2 = st.columns(2)
    with col1:
        linker_smiles = st.text_input("Linker SMILES")
    with col2:
        glycan_smiles = st.text_input("Glycan SMILES")
    bond_idx = st.number_input("Bonding Residue Index", value=50)
    
    if st.button("Save Antigen Info"):
        st.session_state.last_antigen_prot = prot_seq
        st.session_state.last_linker_smiles = linker_smiles
        st.session_state.last_smiles = glycan_smiles
        st.session_state.last_bond_idx = bond_idx
        st.success("抗原情報を保存しました。")

with tab5:
    st.header("🎨 Antibody Candidate Ranking")
    designer = AntibodyDesigner()
    motif = st.selectbox("Target Motif", ["Tn Antigen"])
    ranked = designer.get_ranked_candidates(motif)
    st.table(pd.DataFrame(ranked)[["name", "Score"]])
    
    selected = st.selectbox("Select Candidate to Export", [r["name"] for r in ranked])
    best = next(r for r in ranked if r["name"] == selected)
    
    if st.button("Generate AF3 Full JSON"):
        if st.session_state.last_antigen_prot:
            wf = GlycoConjugateWorkflow(job_name)
            full_json = wf.create_full_complex_json(
                job_name, st.session_state.last_antigen_prot, 
                st.session_state.last_smiles, 
                st.session_state.get('last_linker_smiles', ""),
                st.session_state.last_bond_idx, 
                best["H_AA"], best["L_AA"]
            )
            st.download_button("Download Full JSON", json.dumps(full_json, indent=2), f"{job_name}_full.json", "application/json")
        else:
            st.warning("先に「🧬 抗原デザイン」タブで情報を保存してください。")