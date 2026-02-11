import streamlit as st
import pandas as pd
import json
from your_module import GlycoConjugateWorkflow, ComplexBuilder, AntibodyGraftingEngine, CDRPredictor, HotSpotAnalyzer

st.set_page_config(page_title="GlycoVaccine Studio v2.0", layout="wide", page_icon="🧪")
st.title("🧪 GlycoVaccine Studio v2.0")

# --- Session State ---
if 'ant_seq' not in st.session_state: st.session_state.ant_seq = ""
if 'l_smi' not in st.session_state: st.session_state.l_smi = ""
if 'g_smi' not in st.session_state: st.session_state.g_smi = ""
if 'h_seq' not in st.session_state: st.session_state.h_seq = ""
if 'l_seq' not in st.session_state: st.session_state.l_seq = ""

with st.sidebar:
    st.header("⚙️ Settings")
    platform = st.radio("Target Platform", ["AlphaFold Server (Web)", "AlphaFold 3 (Standalone)"])
    job_name = st.text_input("Job Name", "TACA_Project")

tab1, tab2, tab3, tab4 = st.tabs(["🧬 1. 複合体作製", "🎨 2. 抗体設計", "🛡️ 3. 解析 & JSON出力", "🔥 4. Hot Spot"])

with tab1:
    st.header("🧬 Antigen-Glycan Complex Builder")
    prot = st.text_area("Carrier Protein Sequence", value=st.session_state.ant_seq)
    col1, col2 = st.columns(2)
    with col1: l_smi = st.text_input("Linker SMILES", value=st.session_state.l_smi)
    with col2: g_smi = st.text_input("Glycan SMILES", value=st.session_state.g_smi)
    if st.button("Build Complex"):
        st.session_state.ant_seq, st.session_state.l_smi, st.session_state.g_smi = prot, l_smi, g_smi
        builder = ComplexBuilder()
        pdb_data = builder.build_antigen_pdb(prot, l_smi, g_smi)
        st.download_button("Download PDB", pdb_data, "antigen_complex.pdb")

with tab2:
    st.header("🎨 Antibody Engineering (Trastuzumab Grafting)")
    if st.button("Graft CDRs onto Trastuzumab"):
        pred = CDRPredictor()
        h_cdrs, l_cdrs = pred.predict(st.session_state.g_smi)
        engine = AntibodyGraftingEngine()
        st.session_state.h_seq, st.session_state.l_seq = engine.graft(h_cdrs, l_cdrs)
        st.success("CDR を移植しました。")
        st.text_area("Heavy Chain", st.session_state.h_seq)
        st.text_area("Light Chain", st.session_state.l_seq)

with tab3:
    st.header("🛡️ Output for AF3 & CueMol2")
    if st.session_state.ant_seq and st.session_state.h_seq:
        col3, col4 = st.columns(2)
        with col3:
            if st.button("Generate AF3 JSON"):
                wf = GlycoConjugateWorkflow(job_name)
                full_json = wf.create_full_complex_json(
                    job_name, st.session_state.ant_seq, st.session_state.g_smi, 
                    st.session_state.l_smi, 50, st.session_state.h_seq, st.session_state.l_seq, 
                    mode=platform
                )
                st.download_button("Download JSON", json.dumps(full_json, indent=2), f"{job_name}_full.json")
        with col4:
            if st.button("Generate PDB for CueMol2"):
                                builder = ComplexBuilder()
            merged = builder.merge_for_cuemol(st.session_state.ant_seq, st.session_state.h_seq, st.session_state.l_seq)
            st.download_button("Download Merged PDB", merged, "complex_for_cuemol.pdb")
    else:
        st.warning("先に Tab 1 と 2 でデータを作成してください。")

with tab4:
    st.header("🔥 Hot Spot Analysis")
    uploaded = st.file_uploader("Upload Antibody Model")
    if uploaded:
        ana = HotSpotAnalyzer()
        df = ana.analyze(uploaded.name)
        st.table(df)
        st.bar_chart(df.set_index("Residue"))