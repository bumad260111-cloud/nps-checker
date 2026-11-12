import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import MolToImage

st.set_page_config(page_title="Экспертиза НПВ — Постановление № 1213", layout="wide")

st.title("🔬 Автоматизация отнесения веществ к производным")
st.caption("Проверка по структурно-групповым признакам и автоопределение наименований (ПКМУ № 770 / № 1213)")

# Боковая панель для ввода данных
st.sidebar.header("Ввод структуры")
smiles_input = st.sidebar.text_input("Введите SMILES соединения:", "O=C(c1ccccc1)C(NC)CC")
search_btn = st.sidebar.button("Проверить вещество")

# Базовые SMARTS-шаблоны
PATTERNS = {
    "Производные синтетических катинонов": "c1ccccc1CC(N)C=O",
    "Производные фенетиламинов": "c1ccccc1CCN",
    "Синтетические канабиноиды (Индолы/Индазолы)": "c1ccc2c(c1)c(c[nH]2)",
    "Производные фентанилов": "N1CCC(CC1)(N(C(=O)CC)c2ccccc2)c3ccccc3",
    "Производные триптаминов": "c1ccc2c(c1)c(c[nH]2)CCN"
}

def get_pubchem_info(smiles):
    """Запрос официального названия и CAS из PubChem по SMILES"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName,Title/json"
        response = requests.get(url, timeout=3)
        if response.status_code == 200:
            data = response.json()
            props = data['PropertyTable']['Properties'][0]
            iupac_name = props.get('IUPACName', 'Не найдено')
            title = props.get('Title', 'Не найдено')
            return title, iupac_name
    except Exception:
        pass
    return "Сведения уточняются", "Не удалось получить данные из сети"

def analyze_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, "Ошибка: неверный формат SMILES. Проверьте правильность строки."

    mw = Descriptors.MolWt(mol)
    heavy_atoms_count = mol.GetNumHeavyAtoms()
    
    matched_groups = []
    for group_name, smarts in PATTERNS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            matched_groups.append(group_name)

    # Получаем информацию из открытых баз данных
    common_name, iupac_name = get_pubchem_info(smiles)

    return {
        "mol": mol,
        "mw": round(mw, 2),
        "heavy_atoms": heavy_atoms_count,
        "groups": matched_groups,
        "common_name": common_name,
        "iupac_name": iupac_name
    }, None

if smiles_input:
    res, err = analyze_molecule(smiles_input)
    if err:
        st.error(err)
    else:
        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.subheader("Визуализация структуры")
            img = MolToImage(res["mol"], size=(300, 300))
            st.image(img, use_column_width=True)
            st.metric("Молекулярная масса (а.о.м.)", res["mw"])
            st.metric("Количество тяжелых атомов", res["heavy_atoms"])

        with col2:
            st.subheader("Результаты анализа и идентификация")
            st.info(, icon="ℹ️")
            st.write(f"**Наименование / Синоним:** {res['common_name']}")
            st.write(f"**Систематическое название (IUPAC):** {res['iupac_name']}")
            
            st.markdown("---")
            if res["groups"]:
                for g in res["groups"]:
                    st.success(f"Обнаружено структурное соответствие группе: **{g}**")
                
                st.markdown("---")
                st.subheader("📝 Проект вывода эксперта:")
                st.code(
                    f"Исследуемое вещество ({res['common_name']}) содержит базовое ядро и структурные "
                    f"признаки, характерные для группы «{res['groups'][0]}».\n"
                    f"Согласно Постановлению КМУ № 1213 от 02.10.2026 о внесении изменений в ПКМУ № 770, "
                    f"данное соединение подпадает под критерии отнесения к производным подконтрольных веществ.",
                    language="text"
                )
            else:
                st.warning("Прямых совпадений с базовыми ядрами подконтрольных групп в рамках текущих шаблонов не обнаружено. Требуется ручной анализ.")
