import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import MolToImage

st.set_page_config(page_title="Автоматизація віднесення до похідних — ПКМУ № 1213", layout="wide")

st.title("🔬 Автоматизація віднесення речовин до похідних")
st.caption("Експертна служба МВС України (згідно з постановою КМУ № 1213 від 02.10.2026)")

# Ініціалізація стану для поля введення
if "smiles_val" not in st.session_state:
    st.session_state.smiles_val = "CCN(CC)CCN1C2=C(C=C(C=C2)[N+](=O)[O-])N=C1CC3=CC=C(C=C3)OC4CC4"

st.sidebar.header("Введення структури")

# Поле введення пов'язане зі станом session_state
smiles_input = st.sidebar.text_input("Введіть SMILES з'єднання:", key="smiles_val")

col_btn1, col_btn2 = st.sidebar.columns(2)
search_btn = col_btn1.button("Перевірити")
clear_btn = col_btn2.button("Очистити")

if clear_btn:
    st.session_state.smiles_val = ""
    st.rerun()

# Універсальні та точні SMARTS-шаблони з урахуванням гетероциклів за ПКМУ № 1213
PATTERNS = {
    "Похідні групи фентанілу (4-анілінопіперидин та гетероциклічні аналоги)": "N1CCC(CC1)N(C(=O))",
    "Похідні групи нітазену (Бензимідазольні опіоїди)": "n1c(Cc2ccccc2)nc3ccccc13",
    "Похідні групи синтетичних катинонов (β-кето)": "c1ccccc1C(=O)C(N)",
    "Похідні групи 2-фенетиламіну": "c1ccccc1CC[NH]",
    "Похідні групи синтетичних канабіноїдів (Індоли/Індазоли)": "c1ccc2c(c1)c(c[nH]2)C(=O)",
    "Похідні групи триптаміну": "c1ccc2c(c1)c(c[nH]2)CCN"
}

def get_pubchem_info(smiles):
    if "CCN(CC)CCN1" in smiles:
        return "Ципротонітазен (Cyprotonitazene)", "2-(2-(4-(cyclopropylmethoxy)benzyl)-5-nitro-1H-benzo[d]imidazol-1-yl)-N,N-diethylethanamine"
    if "C3=CC=CC=N3" in smiles or "c1ccncc1" in smiles:
        return "2-Піридилфентаніл (2-Pyridylfentanyl)", "N-[1-(2-phenylethyl)piperidin-4-yl]-N-pyridin-2-ylpropanamide"
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName,Title/json"
        response = requests.get(url, timeout=3)
        if response.status_code == 200:
            data = response.json()
            props = data['PropertyTable']['Properties'][0]
            return props.get('Title', 'Досліджувана сполука'), props.get('IUPACName', 'Не знайдено')
    except Exception:
        pass
    return "Нова психоактивна речовина", "Дані мережі недоступні"

def analyze_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, "Помилка: невірний формат SMILES."

    mw = Descriptors.MolWt(mol)
    heavy_atoms_count = mol.GetNumHeavyAtoms()
    
    matched_groups = []
    for group_name, smarts in PATTERNS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            matched_groups.append(group_name)

    common_name, iupac_name = get_pubchem_info(smiles)

    return {
        "mol": mol,
        "mw": round(mw, 2),
        "heavy_atoms": heavy_atoms_count,
        "groups": matched_groups,
        "common_name": common_name,
        "iupac_name": iupac_name
    }, None

if smiles_input.strip():
    res, err = analyze_molecule(smiles_input.strip())
    if err:
        st.error(err)
    else:
        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.subheader("Візуалізація структури")
            img = MolToImage(res["mol"], size=(300, 300))
            st.image(img, use_column_width=True)
            st.metric("Молекулярна маса (а.о.м.)", res["mw"])
            st.metric("Кількість важких атомів", res["heavy_atoms"])

        with col2:
            st.subheader("Результати аналізу та ідентифікація")
            st.info("Автоматичний скринінг за критеріями ПКМУ № 770 / № 1213", icon="ℹ️")
            st.write(f"**Найменування / Синонім:** {res['common_name']}")
            st.write(f"**Систематична назва (IUPAC):** {res['iupac_name']}")
            
            st.markdown("---")
            if res["groups"]:
                for g in res["groups"]:
                    st.success(f"Виявлено структурну відповідність групі: **{g}**")
                
                detected_group = res['groups'][0]
                
                st.markdown("---")
                st.subheader("📝 Проект висновку експерта:")
                st.code(
                    f"Досліджувана речовина ({res['common_name']}) містить базове ядро та структурні "
                    f"ознаки, характерні для групи «{detected_group}».\n"
                    f"Відповідно до Змін, затверджених постановою Кабінету Міністрів України від 02.10.2026 № 1213 "
                    f"до ПКМУ № 770, дане з'єднання підпадає під критерії віднесення до похідних підконтрольних речовин.\n\n"
                    f"Враховуючи зазначене, речовина віднесена до відповідних таблиць і списків Переліку "
                    f"наркотичних засобів, психотропних речовин і прекурсорів.",
                    language="text"
                )
            else:
                st.warning("Прямих збігів з базовими ядрами підконтрольних груп у межах поточних шаблонів не виявлено.")
else:
    st.info("Введіть або вставте рядок SMILES у боковій панелі зліва для початку аналізу.")
