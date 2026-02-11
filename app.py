import streamlit as st
import pandas as pd
from your_module import AntibodyGraftingEngine, ComplexBuilder, CDRPredictor, HotSpotAnalyzer

st.set_page_config(page_title="GlycoVaccine Studio v2.0", layout="wide")
st.title("🧪 GlycoVaccine Studio v2.0")

# --- Session State ---
if 'antigen_complex' not in st.session_state: st.session_state.antigen_complex = None
if 'engineered_antibody' not in st.session_state: st.session_state.engineered_antibody = None

tab1, tab2, tab3, tab4 = st.tabs([
    "🧬 1. 複合体作製", "🎨 2. 抗体エンジニアリング", "🛡️ 3. 結合図・CueMol出力", "🔥 4. Hot Spot解析"
])

# --- Tab 1: 複合体作製 ---
with tab1:
    st.header("🧬 Antigen-Glycan Complex Builder")
    prot_seq = st.text_area("Carrier Protein (CRM197等)")
    col1, col2 = st.columns(2)
    with col1: linker_smi = st.text_input("Linker SMILES")
    with col2: glycan_smi = st.text_input("Glycan SMILES")
    
    if st.button("Build & Output Complex"):
        builder = ComplexBuilder()
        complex_data = builder.build_antigen_glycan_cif(prot_seq, linker_smi, glycan_smi)
        st.session_state.antigen_complex = complex_data
        st.download_button("Download Complex (CIF)", complex_data, "antigen_complex.cif")
        st.success("複合体情報を生成しました。")

# --- Tab 2: 抗体エンジニアリング ---
with tab2:
    st.header("🎨 CDR Grafting (Trastuzumab Template)")
    if st.button("Predict CDRs for Target Glycan"):
        predictor = CDRPredictor()
        h_cdrs, l_cdrs = predictor.predict_for_glycan(glycan_smi)
        
        graft_engine = AntibodyGraftingEngine()
        h_full, l_full = graft_engine.graft_cdrs(h_cdrs, l_cdrs)
        
        st.session_state.engineered_antibody = {"H": h_full, "L": l_full}
        st.subheader("Grafted Antibody Sequences")
        st.text_area("Heavy Chain (Trastuzumab FR + Predicted CDRs)", h_full)
        st.text_area("Light Chain (Trastuzumab FR + Predicted CDRs)", l_full)
        st.download_button("Download Antibody Sequence", f">H_chain\n{h_full}\n>L_chain\n{l_full}", "antibody.fasta")

# --- Tab 3: 結合図 & CueMol2出力 ---
with tab3:
    st.header("🛡️ Visualizer & CueMol2 Export")
    st.info("作成した複合体と抗体を統合して出力します。")
    if st.session_state.antigen_complex and st.session_state.engineered_antibody:
        # ここで3D可視化コンポーネントを表示（前回同様の3Dmol.jsを利用）
        st.button("Combine for CueMol2")
        st.download_button("Download for CueMol2 (PDB)", "PDB_DATA_HERE", "complex_for_cuemol.pdb")
    else:
        st.warning("先にTab 1とTab 2でデータを作成してください。")

# --- Tab 4: Hot Spot解析 ---
with tab4:
    st.header("🔥 Hot Spot Analysis")
    uploaded_ab = st.file_uploader("Upload Antibody Model")
    if uploaded_ab:
        analyzer = HotSpotAnalyzer()
        df = analyzer.analyze(uploaded_ab.name)
        st.table(df)
        st.bar_chart(df.set_index("Residue"))