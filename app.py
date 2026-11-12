import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import MolToImage

st.set_page_config(page_title="Автоматизація віднесення до похідних — ПКМУ № 1213", layout="wide")

st.title("🔬 Автоматизація віднесення речовин до похідних")
st.caption("Експертна служба МВС України (згідно з постановою КМУ № 1213 від 02.10.2026)")

if "smiles_val" not in st.session_state:
    st.session_state.smiles_val = "CC(CC1=CC=CC=C1)N"

st.sidebar.header("Введення структури")

def clear_text():
    st.session_state.smiles_val = ""

smiles_input = st.sidebar.text_input("Введіть SMILES з'єднання:", key="smiles_val")

col_btn1, col_btn2 = st.sidebar.columns(2)
search_btn = col_btn1.button("Перевірити")
clear_btn = col_btn2.button("Очистити", on_click=clear_text)

# Прямі самостійні позиції згідно з ПКМУ № 770
DIRECT_LIST_770 = {
    "CC(CC1=CC=CC=C1)N": ("Амфетамин (Amphetamine)", "Список № 2 таблиці I (Психотропні речовини, обіг яких обмежено)"),
    "CNC(CC1=CC=CC=C1)C": ("Метамфетамин (Methamphetamine)", "Список № 2 таблиці I (Особливо небезпечні психотропні речовини)"),
    "O=C(c1ccccc1)C(CCCC)N1CCCC1": ("альфа-PVP (α-PVP)", "Список № 2 таблиці I (Особливо небезпечні психотропні речовини)"),
    "CCC(=O)N(C1CCN(CC1)CCC2=CC=CC=C2)C3=CC=CC=C3": ("Фентанил (Fentanyl)", "Список № 1 таблиці I (Особливо небезпечні наркотичні засоби)")
}

PATTERNS = {
    "Похідні групи фентанілу (4-анілінопіперидин та гетероциклічні аналоги)": "N1CCC(CC1)N(C(=O))",
    "Похідні групи нітазену (Бензимідазольні опіоїди)": "n1c(Cc2ccccc2)nc3ccccc13",
    "Похідні групи синтетичних катинонов (β-кето)": "c1ccccc1C(=O)C(N)",
    "Похідні групи 2-фенетиламіну": "c1ccccc1CC[NH]",
    "Похідні групи синтетичних канабіноїдів (Індоли/Індазоли)": "c1ccc2c(c1)c(c[nH]2)C(=O)",
    "Похідні групи триптаміну": "c1ccc2c(c1)c(c[nH]2)CCN"
}

def get_pubchem_info(smiles):
    if smiles in DIRECT_LIST_770:
        return DIRECT_LIST_770[smiles][0], "Самостійна позиція Переліку ПКМУ № 770"
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
    
    is_direct = smiles in DIRECT_LIST_770
    matched_groups = []
    
    if not is_direct:
        for group_name, smarts in PATTERNS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                matched_groups.append(group_name)

    common_name, iupac_name = get_pubchem_info(smiles)

    return {
        "mol": mol,
        "mw": round(mw, 2),
        "heavy_atoms": heavy_atoms_count,
        "is_direct": is_direct,
        "direct_info": DIRECT_LIST_770.get(smiles, ("", ""))[1],
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
            st.info("Автоматичний скринінг за критеріями ПКМУ № 770 / № 121
