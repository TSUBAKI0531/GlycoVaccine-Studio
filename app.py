import streamlit as st
import pandas as pd
import json
from your_module import GlycoConjugateWorkflow, AntibodyDesigner

# --- Session State 初期化 ---
if 'last_antigen_prot' not in st.session_state: st.session_state.last_antigen_prot = ""
if 'last_smiles' not in st.session_state: st.session_state.last_smiles = ""
if 'last_linker_smiles' not in st.session_state: st.session_state.last_linker_smiles = ""
if 'last_bond_idx' not in st.session_state: st.session_state.last_bond_idx = 50

st.set_page_config(page_title="GlycoVaccine Studio", layout="wide", page_icon="🧪")
st.title("🧪 GlycoVaccine Studio v1.5")

with st.sidebar:
    st.header("⚙️ Settings")
    platform = st.radio("Target Platform", ["AlphaFold Server (Web)", "AlphaFold 3 (Standalone)"])
    job_name = st.text_input("Job Name", "TACA_Project")

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["🧬 抗原デザイン", "📊 モデル解析", "🛡️ 抗体解析", "🔥 Hot Spot", "🎨 抗体デザイン", "📋 履歴"])

# --- Tab 1: 抗原デザイン ---
with tab1:
    st.header("🧬 Antigen Design with Linker")
    prot_seq = st.text_area("Antigen Protein Sequence (Carrier)", value=st.session_state.last_antigen_prot)
    col1, col2 = st.columns(2)
    with col1:
        linker_smi = st.text_input("Linker SMILES", value=st.session_state.last_linker_smiles)
    with col2:
        glycan_smi = st.text_input("Glycan SMILES", value=st.session_state.last_smiles)
    bond_idx = st.number_input("Bonding Residue Index", value=st.session_state.last_bond_idx)
    
    if st.button("Save Antigen Info"):
        st.session_state.last_antigen_prot = prot_seq
        st.session_state.last_linker_smiles = linker_smi
        st.session_state.last_smiles = glycan_smi
        st.session_state.last_bond_idx = bond_idx
        st.success("抗原情報を保存しました。")

# --- Tab 5: 抗体デザイン（ここでTypeErrorが発生していた箇所） ---
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
            # 引数の数をyour_module.pyの定義(7個+keyword)と完全に一致させる
            full_json = wf.create_full_complex_json(
                job_name, 
                st.session_state.last_antigen_prot, 
                st.session_state.last_smiles, 
                st.session_state.last_linker_smiles, 
                st.session_state.last_bond_idx, 
                best["H_AA"], 
                best["L_AA"], 
                platform=platform
            )
            st.download_button("Download JSON", json.dumps(full_json, indent=2), f"{job_name}_full.json")
        else:
            st.warning("先に「🧬 抗原デザイン」タブで情報を保存してください。")