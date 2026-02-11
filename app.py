import streamlit as st
import streamlit.components.v1 as components
from your_module import ComplexBuilder

# --- 3D 可視化コンポーネント ---
def show_3d_viewer(pdb_text):
    """3Dmol.js を使用して PDB データを描画する"""
    html_code = f"""
    <div id="container" style="height: 400px; width: 100%;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        $(function() {{
            let viewer = $3Dmol.createViewer($("#container"), {{backgroundColor: "white"}});
            viewer.addModel(`{pdb_text}`, "pdb");
            viewer.setStyle({{cartoon: {{color: 'spectrum'}}}});
            viewer.zoomTo();
            viewer.render();
        }});
    </script>
    """
    components.html(html_code, height=420)

st.set_page_config(page_title="GlycoVaccine Studio", layout="wide")
st.title("🧪 GlycoVaccine Studio")

# タブの作成（まずは Tab 1 から）
tab1, = st.tabs(["🧬 1. 複合体作製"])

with tab1:
    st.header("Antigen-Glycan Complex Builder")
    
    # 入力エリア
    prot_seq = st.text_area("Carrier Protein Sequence (例: CRM197 等)", height=150)
    col1, col2 = st.columns(2)
    with col1:
        linker_smi = st.text_input("Linker SMILES", "NCCCO")
    with col2:
        glycan_smi = st.text_input("Glycan SMILES", "NCCCO[C@@H]1[C@@H](NC(C)=O)[C@@H](O)O[C@H](CO)[C@H]1O")
    
    if st.button("複合体を作製して表示"):
        if prot_seq:
            builder = ComplexBuilder()
            # 座標データの生成
            pdb_data = builder.build_complex_pdb(prot_seq, linker_smi, glycan_smi)
            
            st.success("タンパク質-糖鎖複合体モデル（主鎖）を生成しました。")
            
            # アプリ上での 3D 表示
            show_3d_viewer(pdb_data)
            
            # CueMol2 用のファイル出力
            st.download_button(
                label="CueMol2 用 PDB をダウンロード",
                data=pdb_data,
                file_name="antigen_complex.pdb",
                mime="chemical/x-pdb"
            )
        else:
            st.warning("タンパク質配列を入力してください。")