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

# --- Session State ---
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
    prot_seq = st.text_area("Antigen Protein Sequence")
    smiles = st.text_input("Glycan SMILES")
    bond_idx = st.number_input("Bonding Residue Index", value=1)
    if st.button("Save Antigen Info"):
        st.session_state.last_antigen_prot = prot_seq
        st.session_state.last_smiles = smiles
        st.session_state.last_bond_idx = bond_idx
        st.success("抗原情報を保存しました。")

with tab2:
    uploaded = st.file_uploader("Upload CIF models", accept_multiple_files=True)
    if uploaded:
        wf = GlycoConjugateWorkflow(job_name)
        for f in uploaded:
            content = f.read().decode("utf-8")
            with open(f.name, "w") as tmp: tmp.write(content)
            res = wf._calculate_sasa(f.name)
            st.write(f"Model: {f.name}, SASA: {res['glycan_sasa']:.2f}")
            st.session_state.analysis_history.append({"Job": job_name, "Type": "SASA", "Target": f.name, "Score": res['glycan_sasa']})

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

with tab6:
    if st.session_state.analysis_history:
        df_h = pd.DataFrame(st.session_state.analysis_history)
        st.dataframe(df_h)
        output = BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            df_h.to_excel(writer, index=False)
        st.download_button("Download Excel", output.getvalue(), "history.xlsx")