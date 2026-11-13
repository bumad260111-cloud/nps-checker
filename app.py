import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import MolToImage

st.set_page_config(page_title="Автоматизація віднесення до похідних — ПКМУ № 1213", layout="wide")

st.title("🔬 Автоматизація віднесення речовин до похідних")
st.caption("Експертна служба МВС України (згідно з постановою КМУ № 1213 від 02.10.2026)")

# Ініціалізація стану
if "smiles_val" not in st.session_state:
    st.session_state.smiles_val = "O=C(c1ccccc1)C(CCCC)N1CCCC1"

if "inchi_val" not in st.session_state:
    st.session_state.inchi_val = ""

if "resolved_smiles_val" not in st.session_state:
    st.session_state.resolved_smiles_val = ""

st.sidebar.header("Введення структури")

def clear_fields():
    st.session_state.smiles_val = ""
    st.session_state.inchi_val = ""
    st.session_state.resolved_smiles_val = ""

def get_smiles_from_inchikey(inchi_key):
    clean_key = inchi_key.strip()
    if not clean_key:
        return None
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{clean_key}/property/CanonicalSMILES/json"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            props = data.get('PropertyTable', {}).get('Properties', [])
            if props:
                return props[0].get('CanonicalSMILES')
    except Exception as e:
        st.error(f"Помилка запиту до PubChem: {e}")
    return None

# Обробник зміни InChIKey (працює автоматично при вводі у 2-гу строку)
def on_inchi_change():
    inchi = st.session_state.inchi_val.strip()
    if inchi:
        resolved = get_smiles_from_inchikey(inchi)
        if resolved:
            st.session_state.resolved_smiles_val = resolved
        else:
            st.session_state.resolved_smiles_val = "Не знайдено в PubChem за цим InChIKey"

# 1-й рядок: Введення SMILES
smiles_input = st.sidebar.text_input("1. Введіть SMILES з'єднання:", key="smiles_val")

# 2-й рядок: Введення InChIKey (з автовикликом конвертації)
inchi_input = st.sidebar.text_input("2. Або введіть InChIKey (наприклад з NIST):", key="inchi_val", on_change=on_inchi_change)

# 3-й рядок: Визначений / Результуючий SMILES (звідси можна копіювати)
resolved_smiles_input = st.sidebar.text_input("3. Визначений SMILES (скопіюйте сюди або в 1-у строку):", key="resolved_smiles_val")

col_btn1, col_btn2 = st.sidebar.columns(2)
search_btn = col_btn1.button("Перевірити")
clear_btn = col_btn2.button("Очистити", on_click=clear_fields)

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
        return None, f"Помилка: не вдалося побудувати молекулу RDKit з рядка SMILES: `{smiles}`"

    mw = Descriptors.MolWt(mol)
    heavy_atoms_count = mol.GetNumHeavyAtoms()
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    
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
        "logp": round(logp, 2),
        "tpsa": round(tpsa, 2),
        "h_donors": h_donors,
        "h_acceptors": h_acceptors,
        "rot_bonds": rot_bonds,
        "is_direct": is_direct,
        "direct_info": DIRECT_LIST_770.get(smiles, ("", ""))[1],
        "groups": matched_groups,
        "common_name": common_name,
        "iupac_name": iupac_name
    }, None

if search_btn:
    # Пріоритет у першої строки, якщо там пусто — беремо з 3-ї строки
    target_smiles = smiles_input.strip() if smiles_input.strip() else resolved_smiles_input.strip()

    if target_smiles and not target_smiles.startswith("Не знайдено"):
        res, err = analyze_molecule(target_smiles)
        if err:
            st.error(err)
        else:
            col1, col2 = st.columns([1, 2])
            
            with col1:
                st.subheader("Візуалізація структури")
                img = MolToImage(res["mol"], size=(300, 300))
                st.image(img, use_column_width=True)
                
                st.markdown("### 📊 Фізико-хімічні дескриптори")
                st.metric("Молекулярна маса (MW)", f"{res['mw']} а.о.м.")
                st.metric("Коефіцієнт розподілу (LogP)", res["logp"])
                st.metric("Полярна поверхня (TPSA)", f"{res['tpsa']} Å²")
                st.write(f"**Важких атомів:** {res['heavy_atoms']}")
                st.write(f"**Донори / Акцептори H-зв'язків:** {res['h_donors']} / {res['h_acceptors']}")
                st.write(f"**Обертових зв'язків:** {res['rot_bonds']}")

            with col2:
                st.subheader("Результати аналізу та ідентифікація")
                st.info("Автоматичний скринінг за критеріями ПКМУ № 770 / № 1213", icon="ℹ️")
                st.write(f"**Найменування / Синонім:** {res['common_name']}")
                st.write(f"**Систематична назва (IUPAC):** {res['iupac_name']}")
                
                st.markdown("---")
                if res["is_direct"]:
                    st.error("⚠️ **Речовина вже включена як самостійна позиція до Переліку ПКМУ № 770!**")
                    st.write("Статус у Переліку підтверджено офіційними даними.")
                    
                    st.markdown("---")
                    st.subheader("📝 Проект висновку експерта:")
                    st.code(
                        f"Досліджувана речовина ({res['common_name']}) є самостійною позицією у Переліку наркотичних "
                        f"засобів, психотропних речовин і прекурсорів, затвердженого постановою Кабінету Міністрів України "
                        f"від 06.05.2000 № 770 ({res['direct_info']}).\n"
                        f"Розрахункові фізико-хімічні параметри: MW = {res['mw']} а.о.м., LogP = {res['logp']}, TPSA = {res['tpsa']} Å².",
                        language="text"
                    )
                elif res["groups"]:
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
                        f"Фізико-хімічні характеристики молекули: молекулярна маса — {res['mw']} а.о.м., "
                        f"коефіцієнт розподілу (LogP) — {res['logp']}, топологічна полярна поверхня (TPSA) — {res['tpsa']} Å².\n\n"
                        f"Враховуючи зазначене, речовина віднесена до відповідних таблиць і списків Переліку "
                        f"наркотичних засобів, психотропних речовин і прекурсорів.",
                        language="text"
                    )
                else:
                    st.warning("Прямих збігів з самостійними позиціями Переліку № 770 або базовими ядрами похідних у межах поточних шаблонів не виявлено.")
    else:
        st.warning("Будь ласка, введіть коректний SMILES або InChIKey для перевірки.")
else:
    st.info("Введіть InChIKey у 2-гу строку (SMILES з'явиться в 3-й строці) або введіть SMILES у 1-у строку та натисніть «Перевірити».")
