import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import MolToImage

st.set_page_config(page_title="Экспертиза НПВ — Постановление № 1213", layout="wide")

st.title("🔬 Автоматизация отнесения веществ к производным")
st.caption("Определение принадлежности согласно Постановлению КМУ № 1213 от 02.10.2026")

# Боковая панель для ввода данных
st.sidebar.header("Ввод структуры")
smiles_input = st.sidebar.text_input("Введите SMILES соединения:", "CCN(CC)C(=O)C1=CC=CC=C1")
search_btn = st.sidebar.button("Проверить вещество")

# Базовые SMARTS-шаблоны ядер
PATTERNS = {
    "Fentanyls": "[NX3][CX4][CX4][NX3]",  # Упрощенный фрагмент фентанилового ядра
    "Nitazenes": "c1ccc2c(c1)[nH]c(n2)CC3ccc(cc3)",  # Бензимидазольное ядро
    "Cathinones / 2-Phenethylamines": "c1ccccc1CCN",  # Фенетиламиновый скелет
    "Tryptamines": "c1ccc2c(c1)c(c[nH]2)CCN",  # Индольный скелет триптамина
    "Synthetic Cannabinoids": "c1ccc2c(c1)c(c[nH]2)"  # Базовое ядро индола/индазола
}

def analyze_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, "Ошибка: неверный формат SMILES"

    mw = Descriptors.MolWt(mol)
    
    # Подсчет атомов в непрерывной цепи (без водорода)
    heavy_atoms_count = mol.GetNumHeavyAtoms()
    
    matched_groups = []
    for group_name, smarts in PATTERNS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            matched_groups.append(group_name)

    return {
        "mol": mol,
        "mw": round(mw, 2),
        "heavy_atoms": heavy_atoms_count,
        "groups": matched_groups
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
            st.subheader("Результаты анализа")
            if res["groups"]:
                for g in res["groups"]:
                    st.success(f"Обнаружено соответствие группе: **{g}**")
                
                st.markdown("---")
                st.subheader("📝 Проект вывода эксперта:")
                st.code(
                    f"Исследуемое вещество является производным группы {res['groups'][0]}.\n"
                    f"Согласно Постановлению КМУ № 1213 от 02.10.2026, вещество отнесено "
                    f"к спискам подконтрольных веществ (Таблица I, Список № 1/2).",
                    language="text"
                )
            else:
                st.info("Признаков совпадения с базовыми подконтрольными группами не обнаружено.")
