import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
from your_module import StructureMerger, AntibodyDockingWorkflow, AntibodyDesigner

# --- 3D可視化関数 ---
def show_3d_viewer(cif_text):
    html_code = f"""
    <div id="container" style="height: 500px; width: 100%;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        $(function() {{
            let viewer = $3Dmol.createViewer($("#container"), {{backgroundColor: "white"}});
            viewer.addModel(`{cif_text.replace("`", "\\`").replace("$", "\\$")}`, "mcif");
            viewer.setStyle({{cartoon: {{color: 'spectrum'}}}});
            viewer.setStyle({{hetflag: true}}, {{stick: {{radius: 0.3}}}});
            viewer.zoomTo(); viewer.render();
        }});
    </script>"""
    components.html(html_code, height=520)

st.set_page_config(page_title="GlycoVaccine Studio v1.5", layout="wide")
st.title("🧪 GlycoVaccine Studio v1.5")

tab1, tab2, tab3, tab4, tab5 = st.tabs(["🧬 デザイン", "📊 解析", "🛡️ 結合図作成", "🔥 Hot Spot", "🎨 抗体"])

# --- Tab 3: 結合図作成 (AlphaFold以外の代替手法) ---
with tab3:
    st.header("🛡️ 3D Binding Diagram Generator (Non-AF3)")
    st.info("抗原(糖鎖付)と抗体の個別ファイルを統合して結合図を作成します。")
    
    col1, col2 = st.columns(2)
    with col1:
        antigen_file = st.file_uploader("Upload Antigen-Glycan Model (CIF)", key="ant")
    with col2:
        antibody_file = st.file_uploader("Upload Antibody Model (CIF)", key="ab")
        
    if antigen_file and antibody_file:
        if st.button("Generate & Merge 3D Diagram"):
            with open("antigen.cif", "wb") as f: f.write(antigen_file.read())
            with open("antibody.cif", "wb") as f: f.write(antibody_file.read())
            
            merger = StructureMerger()
            merged_path = merger.merge_structures("antigen.cif", "antibody.cif", "merged_complex.cif")
            
            with open(merged_path, "r") as f:
                merged_content = f.read()
            
            show_3d_viewer(merged_content)
            st.download_button("Download Merged Complex", merged_content, "merged_complex.cif")
            
            # パラトープ解析も実行
            adw = AntibodyDockingWorkflow()
            df = adw.analyze_paratope(merged_path, "H", "L")
            st.subheader("Interface Residue Analysis")
            st.table(df)

# --- Tab 5: 抗体候補 ---
with tab5:
    designer = AntibodyDesigner()
    ranked = designer.get_ranked_candidates("Tn Antigen")
    st.table(pd.DataFrame(ranked)[["name", "Score"]])
    st.code(ranked[0]["H_AA"], language="text") # 配列出力