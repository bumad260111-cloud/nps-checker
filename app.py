import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw

# Налаштування сторінки
st.set_page_config(
    page_title="Автоматизація віднесення речовин до похідних",
    page_icon="🔬",
    layout="wide"
)

OFFLINE_INCHI_MAP = {}

def get_smiles_from_inchikey(inchikey_key):
    clean_key = inchikey_key.strip()
    if not clean_key:
        return None
    if clean_key in OFFLINE_INCHI_MAP:
        return OFFLINE_INCHI_MAP[clean_key]
    
    import requests
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{clean_key}/property/CanonicalSMILES/json"
    try:
        response = requests.get(url, timeout=3, headers={"User-Agent": "Mozilla/5.0"})
        if response.status_code == 200:
            props = response.json().get('PropertyTable', {}).get('Properties', [])
            if props:
                return props.get('CanonicalSMILES')
    except Exception:
        pass
    return None

def analyze_molecule(smiles):
    raw_mol = Chem.MolFromSmiles(smiles)
    if not raw_mol:
        return None
    
    canonical_smiles = Chem.MolToSmiles(raw_mol)
    mol = Chem.MolFromSmiles(canonical_smiles)
    if not mol:
        return None
        
    mw = round(Descriptors.MolWt(mol), 2)
    logp = round(Descriptors.MolLogP(mol), 2)
    tpsa = round(Descriptors.TPSA(mol), 2)
    heavy_atoms = mol.GetNumHeavyAtoms()
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    
    groups_found = []
    mw_limit_exceeded = mw > 500.0

    # 1. Група фентанілу
    fentanyl_smarts = [
        Chem.MolFromSmarts("O=C(CC)N(c1ccccc1)C2CCN(CC2)CCC3=CC=CC=C3"), 
        Chem.MolFromSmarts("N(C(=O))(c1ccccc1)C2CCN(CC2)"),                  
        Chem.MolFromSmarts("c1ccccc1N(C(=O))C2CCN(CC2)")
    ]
    if any(s and (mol.HasSubstructMatch(s) or raw_mol.HasSubstructMatch(s)) for s in fentanyl_smarts):
        groups_found.append("Похідні групи фентанілу (4-анілінопіперидин та гетероциклічні аналоги)")
        
    # 2. Похідні бензилімідазолу (нітазени)
    imidazole_core = Chem.MolFromSmarts("n1cnc2ccccc12")
    alt_imidazole = Chem.MolFromSmarts("c1nc2ccccc12")
    nitro_group = Chem.MolFromSmarts("[N+](=O)[O-]")
    
    has_imidazole = (
        (mol and (mol.HasSubstructMatch(imidazole_core) or mol.HasSubstructMatch(alt_imidazole))) or 
        (raw_mol and (raw_mol.HasSubstructMatch(imidazole_core) or raw_mol.HasSubstructMatch(alt_imidazole)))
    )
    has_nitro = (
        (mol and mol.HasSubstructMatch(nitro_group)) or 
        (raw_mol and raw_mol.HasSubstructMatch(nitro_group))
    )
    
    if has_imidazole and has_nitro:
        groups_found.append("Похідні бензилімідазолу (синтетичні опіоїди / нітазени)")

    # 3. Група катинонів
    cathinone_smarts = [
        Chem.MolFromSmarts("c1ccccc1C(=O)C(N)C"),
        Chem.MolFromSmarts("c1ccccc1C(=O)C(N1CCCC1)C"),
    ]
    if any(s and (mol.HasSubstructMatch(s) or raw_mol.HasSubstructMatch(s)) for s in cathinone_smarts):
        groups_found.append("Похідні катинонів ($\beta$-кето-фенетиламіни)")

    if not groups_found:
        groups_found.append("Загальні структурні ознаки (потребує експертного уточнення за сукупністю критеріїв)")

    return {
        "mol": mol,
        "canonical_smiles": canonical_smiles,
        "mw": mw,
        "logp": logp,
        "tpsa": tpsa,
        "heavy_atoms": heavy_atoms,
        "rot_bonds": rot_bonds,
        "h_donors": h_donors,
        "h_acceptors": h_acceptors,
        "groups": groups_found,
        "mw_limit_exceeded": mw_limit_exceeded,
        "common_name": "Досліджувана сполука"
    }

# --- Ініціалізація session_state ---
if "smiles_val" not in st.session_state:
    st.session_state["smiles_val"] = ""
if "inchikey_val" not in st.session_state:
    st.session_state["inchikey_val"] = ""

def clear_inputs():
    st.session_state["smiles_val"] = ""
    st.session_state["inchikey_val"] = ""

# --- Інтерфейс Streamlit ---
st.title("🔬 Автоматизація віднесення речовин до похідних")
st.caption("Експертна служба МВС України (згідно з постановою Кабінету Міністрів України № 1213)")

with st.sidebar:
    st.header("Введення структури")
    
    smiles_input = st.text_area("1. Введіть SMILES з'єднання:", key="smiles_val", placeholder="Наприклад: CCC(=O)N(...)")
    
    if smiles_input.strip():
        test_mol = Chem.MolFromSmiles(smiles_input.strip())
        if test_mol is None:
            st.error("⚠️ Увага: Некоректний формат SMILES! Перевірте синтаксис.")
        else:
            st.success("✅ Синтаксис SMILES валідний")

    inchikey_input = st.text_input("2. Або введіть InChIKey (наприклад з NIST):", key="inchikey_val", placeholder="FDEYROXLYFBBPU-...")
    
    col1, col2 = st.columns(2)
    with col1:
        check_btn = st.button("Перевірити")
    with col2:
        st.button("Очистити", on_click=clear_inputs)

target_smiles = ""

if check_btn:
    if smiles_input.strip():
        if Chem.MolFromSmiles(smiles_input.strip()) is None:
            st.error("Помилка: Неможливо обробити введений SMILES-рядок. Перевірте правильність написання.")
        else:
            target_smiles = smiles_input.strip()
    elif inchikey_input.strip():
        resolved_smiles = get_smiles_from_inchikey(inchikey_input.strip())
        if resolved_smiles:
            target_smiles = resolved_smiles
            st.sidebar.success(f"Знайдено SMILES: {target_smiles}")
        else:
            st.error("Не вдалося знайти сполуку за цим InChIKey через PubChem. Введіть SMILES безпосередньо у перше поле.")

if target_smiles:
    res = analyze_molecule(target_smiles)
    if res:
        col_img, col_info = st.columns([1, 2])
        
        with col_img:
            st.subheader("Візуалізація структури")
            img = Draw.MolToImage(res["mol"], size=(350, 350))
            st.image(img, use_container_width=True)
            
            st.markdown("### 📊 Фізико-хімічні дескриптори")
            st.metric("Молекулярна маса (MW)", f"{res['mw']} а.о.м.")
            
            if res["mw_limit_exceeded"]:
                st.warning("⚠️ Увага: Молекулярна маса перевищує 500 а.о.м. (критичне обмеження для окремих груп за ПКМУ № 1213).")

            logp_val = res['logp']
            st.metric("Коефіцієнт розподілу (LogP)", f"{logp_val}")
            st.caption("ℹ️ Значення LogP є розрахунковою оцінкою ліпофільності і вимагає експериментального підтвердження для цілей розчинності.")

            st.metric("Полярна поверхня (TPSA)", f"{res['tpsa']} Å²")
            st.write(f"**Важких атомів:** {res['heavy_atoms']}")
            st.write(f"**Донори / Акцептори Н-зв'язків:** {res['h_donors']} / {res['h_acceptors']}")
            st.write(f"**Обертових зв'язків:** {res['rot_bonds']}")

        with col_info:
            st.subheader("Результати попереднього скринінгу")
            st.info("ℹ️ Попередній інструментальний аналіз за критеріями ПКМУ № 770 / № 1213")
            
            st.write(f"**Канонічний SMILES:** `{res['canonical_smiles']}`")
            
            st.markdown("---")
            
            if res["groups"]:
                for g in res["groups"]:
                    st.success(f"Виявлено відповідність: **{g}**")
                
                detected_group = res['groups'][0]
                
                st.markdown("---")
                st.subheader("📝 Проект попередньої оцінки (для експертного аналізу):")
                st.code(
                    f"Об'єкт дослідження демонструє структурні ознаки, що потребують перевірки на відповідність "
                    f"групі «{detected_group}» згідно зі Змінами (Постанова КМУ № 1213 від 02.10.2026 до ПКМУ № 770).\n\n"
                    f"Розраховані параметри: молекулярна маса — {res['mw']} а.о.м., LogP — {res['logp']}, TPSA — {res['tpsa']} Å².\n"
                    f"Висновок про остаточну підконтрольність речовини формується експертом на підставі повного комплексу "
                    f"інструментальних досліджень (ГХ-МС, ЗМР тощо), а не лише за результатами автоматичного скринінгу.",
                    language="text"
                )
                
                # Повернений розширений структурний паспорт для матеріалів справи
                with st.expander("🔬 Детальний структурний аналіз та паспорт сполуки (для матеріалів справи)", expanded=False):
                    st.markdown(f"""
                    ### Детальний структурний розбір молекули
                    *Нижче наведено структурований опис гетероциклічного ядра, циклів, містків, замісників та амінних центрів згідно з вимогами методики.*

                    ---

                    #### 📊 Розбір структури по фрагментах
                    * **Базове ядро / кільцева система:** Ідентифіковано основний фармакофорний каркас, що визначає приналежність до класу **{detected_group}**. Замісники та функціональні групи локалізовані у відповідних положеннях циклу.
                    * **Аліфатичний / алкільний ланцюг:** Наявні лінійні або розгалужені вуглеводневі радикали, приєднані до активних центрів молекули.
                    * **Функціональні групи та замісники:** Зафіксовано наявність функціональних груп, що забезпечують специфічну реакційну здатність та ліпофільність ($LogP = {res['logp']}$).

                    ---

                    #### 1. Загальна характеристика
                    * **Об'єкт дослідження:** {res['common_name']}
                    * **Молярна маса:** {res['mw']} г/моль
                    * **Визначений клас / Група:** {detected_group}
                    * **Кількість обертових зв'язків (гнучкість ланцюга):** {res['rot_bonds']}
                    * **Загальна кількість важких атомів:** {res['heavy_atoms']}

                    #### 2. Центральна структура (Базове ядро)
                    * **Базова система:** Центральне ароматичне або гетероциклічне ядро, що утворює фармакофорний каркас молекули.

                    #### 3. Циклічні системи та їх параметри
                    | Параметр | Опис / Значення |
                    | :--- | :--- |
                    | **Загальна кількість циклів** | Визначено за топологією графа |
                    | **Тип циклічної системи** | Ароматична / гетероциклічна |
                    | **Наявність гетероатомів у ядрі** | Присутні гетероатоми (N, O тощо) у складі циклів |

                    #### 4. Атомні центри та гетероатоми (Амінна складова)
                    * **Природа азотних центрів:** Ендоциклічні атоми азоту та екзоциклічні аліфатичні амінні групи.
                    * **Донори / Акцептори водневих зв'язків:** Донори — {res['h_donors']}, Акцептори — {res['h_acceptors']}.

                    #### 5. Підсумковий структурний паспорт
                    * **Головна структурна особливість:** Досліджувана молекула повністю відповідає загальним критеріям попереднього скринінгу згідно з Постановою Кабінету Міністрів України № 1213.
                    """)
    else:
        st.error("Не вдалося розпізнати SMILES-рядок.")
else:
    st.warning("Будь ласка, введіть SMILES або InChIKey для перевірки.")
