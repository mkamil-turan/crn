"""Protein ve tRNA ölçümlerinin sayısal tabloları.

Her tablo kaynağıyla birlikte verilir; ölçüm fonksiyonlarına parametre olarak
girer ve değiştirilebilir. Doğrulanamamış değerler için
``doc/tables.md`` dosyasına bakın.
"""

# Amino asit grupları. Alifatik, aromatik, küçük, bazik ve asidik gruplar
# yaygın tanımlardır; hidrofobik ve polar gruplar EMBOSS pepstats'ın
# "Non-polar" ve "Polar" sınıflarıdır.
AMINO_ACID_GROUPS = {
    "hydrophobic": "ACFGILMPVWY",
    "polar":       "DEHKNQRST",
    "aliphatic":   "AILV",
    "aromatic":    "FWY",
    "small":       "ACDGNPST",
    "basic":       "HKR",
    "acidic":      "DE",
}

# Kalıntı (amino asit − H2O) kütleleri, Da. Ortalama kütleler doğal izotop
# bolluğuyla, monoizotopik kütleler her elementin en bol izotopuyla
# hesaplanmıştır (ExPASy Compute pI/Mw ve PeptideMass değerleri).
RESIDUE_AVERAGE_MASS = {
    "A": 71.0788, "C": 103.1388, "D": 115.0886, "E": 129.1155, "F": 147.1766,
    "G": 57.0519, "H": 137.1411, "I": 113.1594, "K": 128.1741, "L": 113.1594,
    "M": 131.1926, "N": 114.1038, "P": 97.1167, "Q": 128.1307, "R": 156.1875,
    "S": 87.0782, "T": 101.1051, "V": 99.1326, "W": 186.2132, "Y": 163.1760,
}
RESIDUE_MONOISOTOPIC_MASS = {
    "A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259, "F": 147.06841,
    "G": 57.02146, "H": 137.05891, "I": 113.08406, "K": 128.09496, "L": 113.08406,
    "M": 131.04049, "N": 114.04293, "P": 97.05276, "Q": 128.05858, "R": 156.10111,
    "S": 87.03203, "T": 101.04768, "V": 99.06841, "W": 186.07931, "Y": 163.06333,
}
# Zincirin iki ucundaki H ve OH, yani bir su molekülü, Da.
WATER_AVERAGE_MASS = 18.01524
WATER_MONOISOTOPIC_MASS = 18.01056

# Protein formüllerindeki elementler ve formül tuple'larındaki sıraları.
ELEMENTS = ("C", "H", "N", "O", "S")
# Serbest amino asitlerin formülleri, ELEMENTS sırasıyla (C, H, N, O, S).
AMINO_ACID_FORMULAS = {
    "A": (3, 7, 1, 2, 0),   "C": (3, 7, 1, 2, 1),   "D": (4, 7, 1, 4, 0),
    "E": (5, 9, 1, 4, 0),   "F": (9, 11, 1, 2, 0),  "G": (2, 5, 1, 2, 0),
    "H": (6, 9, 3, 2, 0),   "I": (6, 13, 1, 2, 0),  "K": (6, 14, 2, 2, 0),
    "L": (6, 13, 1, 2, 0),  "M": (5, 11, 1, 2, 1),  "N": (4, 8, 2, 3, 0),
    "P": (5, 9, 1, 2, 0),   "Q": (5, 10, 2, 3, 0),  "R": (6, 14, 4, 2, 0),
    "S": (3, 7, 1, 3, 0),   "T": (4, 9, 1, 3, 0),   "V": (5, 11, 1, 2, 0),
    "W": (11, 12, 2, 2, 0), "Y": (9, 11, 1, 3, 0),
}
# Her peptit bağında ayrılan su, ELEMENTS sırasıyla.
WATER_FORMULA = (0, 2, 0, 1, 0)

# pH 7'de yüklü sayılan kalıntılar (His ihmal edilir).
POSITIVE_RESIDUES = "KR"
NEGATIVE_RESIDUES = "DE"
# Q(pH) hesabında iyonlaşan yan zincirler; uçlar ayrıca "N_term" ve "C_term"
# anahtarlarıyla eklenir.
IONIZABLE_POSITIVE = "HKR"
IONIZABLE_NEGATIVE = "CDEY"

# pKa kümeleri. "N_term" ve "C_term" serbest uç grupların pKa'sıdır;
# "N_term_residue" ve "C_term_residue" uçtaki kalıntıya göre bu değerin
# yerine geçen pKa'lardır (yalnızca Bjellqvist; ExPASy Compute pI/Mw).
PKA_SETS = {
    "EMBOSS": {
        "N_term": 8.6, "C_term": 3.6,
        "K": 10.8, "R": 12.5, "H": 6.5, "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1,
        "N_term_residue": {}, "C_term_residue": {},
    },
    "Lehninger": {
        "N_term": 9.69, "C_term": 2.34,
        "K": 10.5, "R": 12.4, "H": 6.0, "D": 3.86, "E": 4.25, "C": 8.33, "Y": 10.0,
        "N_term_residue": {}, "C_term_residue": {},
    },
    "Bjellqvist": {
        "N_term": 7.5, "C_term": 3.55,
        "K": 10.0, "R": 12.0, "H": 5.98, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0,
        "N_term_residue": {"A": 7.59, "M": 7.0, "S": 6.93, "P": 8.36,
                           "T": 6.82, "V": 7.44, "E": 7.7},
        "C_term_residue": {"D": 4.55, "E": 4.75},
    },
}

# Hidropati / hidrofobiklik ölçekleri (ExPASy ProtScale değerleri).
# Kyte–Doolittle, Hopp–Woods, Eisenberg (normalize konsensüs) ve Janin'de
# büyük değer hidrofobik; Hopp–Woods'ta büyük değer hidrofilik. Wimley–White
# arayüz ölçeğidir (ΔG su -> arayüz, kcal/mol; negatif değer arayüzü tercih
# eder); D, E, K, R yüklü, His nötr değerleriyle.
HYDROPATHY_SCALES = {
    "Kyte-Doolittle": {
        "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2,
        "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5,
        "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
    },
    "Hopp-Woods": {
        "A": -0.5, "C": -1.0, "D": 3.0, "E": 3.0, "F": -2.5, "G": 0.0, "H": -0.5,
        "I": -1.8, "K": 3.0, "L": -1.8, "M": -1.3, "N": 0.2, "P": 0.0, "Q": 0.2,
        "R": 3.0, "S": 0.3, "T": -0.4, "V": -1.5, "W": -3.4, "Y": -2.3,
    },
    "Eisenberg": {
        "A": 0.62, "C": 0.29, "D": -0.90, "E": -0.74, "F": 1.19, "G": 0.48, "H": -0.40,
        "I": 1.38, "K": -1.50, "L": 1.06, "M": 0.64, "N": -0.78, "P": 0.12, "Q": -0.85,
        "R": -2.53, "S": -0.18, "T": -0.05, "V": 1.08, "W": 0.81, "Y": 0.26,
    },
    "Wimley-White": {
        "A": 0.17, "C": -0.24, "D": 1.23, "E": 2.02, "F": -1.13, "G": 0.01, "H": 0.17,
        "I": -0.31, "K": 0.99, "L": -0.56, "M": -0.23, "N": 0.42, "P": 0.45, "Q": 0.58,
        "R": 0.81, "S": 0.13, "T": 0.14, "V": 0.07, "W": -1.85, "Y": -0.94,
    },
    "Janin": {
        "A": 0.3, "C": 0.9, "D": -0.6, "E": -0.7, "F": 0.5, "G": 0.3, "H": -0.1,
        "I": 0.7, "K": -1.8, "L": 0.5, "M": 0.4, "N": -0.5, "P": -0.3, "Q": -0.7,
        "R": -1.4, "S": -0.1, "T": -0.2, "V": 0.6, "W": 0.3, "Y": -0.4,
    },
}
# Boman indeksi için yan zincirlerin sikloheksandan suya transfer serbest
# enerjileri (kcal/mol; Boman 2003, Radzicka–Wolfenden). Hidrofilik kalıntılar
# pozitiftir; Pro için değer yoktur ve 0 alınır.
BOMAN_SCALE = {
    "A": -1.81, "C": -1.28, "D": 8.72, "E": 6.81, "F": -2.98, "G": -0.94, "H": 4.66,
    "I": -4.92, "K": 5.55, "L": -4.92, "M": -2.35, "N": 6.64, "P": 0.0, "Q": 5.54,
    "R": 14.92, "S": 3.40, "T": 2.57, "V": -4.04, "W": -2.33, "Y": 0.14,
}

# 280 nm molar ekstinksiyon katsayıları (M^-1 cm^-1; Pace ve ark. 1995).
# "cystine" bir disülfit bağının (iki Cys) katkısıdır.
EXTINCTION_COEFFICIENTS_280 = {"W": 5500, "Y": 1490, "cystine": 125}
# Kalıntıların kısmi spesifik hacimleri (cm³/g; Cohn ve Edsall 1943).
RESIDUE_SPECIFIC_VOLUMES = {
    "A": 0.74, "C": 0.63, "D": 0.60, "E": 0.66, "F": 0.77, "G": 0.64, "H": 0.67,
    "I": 0.90, "K": 0.82, "L": 0.90, "M": 0.75, "N": 0.62, "P": 0.76, "Q": 0.67,
    "R": 0.70, "S": 0.63, "T": 0.70, "V": 0.86, "W": 0.74, "Y": 0.71,
}
# Kalıntı hacimleri (Å³; Zamyatnin 1972).
RESIDUE_VOLUMES = {
    "A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9, "G": 60.1, "H": 153.2,
    "I": 166.7, "K": 168.6, "L": 166.7, "M": 162.9, "N": 114.1, "P": 112.7, "Q": 143.8,
    "R": 173.4, "S": 89.0, "T": 116.1, "V": 140.0, "W": 227.8, "Y": 193.6,
}
# Gly-X-Gly'deki en büyük çözücüye erişilebilir yüzey alanları (Å²; Tien ve
# ark. 2013, kuramsal değerler).
RESIDUE_MAX_ASA = {
    "A": 129, "C": 167, "D": 193, "E": 223, "F": 240, "G": 104, "H": 224,
    "I": 197, "K": 236, "L": 201, "M": 224, "N": 195, "P": 159, "Q": 225,
    "R": 274, "S": 155, "T": 172, "V": 174, "W": 285, "Y": 263,
}

# Termofilik proteinlerde zenginleşen kalıntılar (Zeldovich ve ark. 2007).
IVYWREL_RESIDUES = "IVYWREL"
# CvP yanlılığı: yüklü (D, E, K, R) ve polar yüksüz (N, Q, S, T) kalıntılar.
CVP_CHARGED = "DEKR"
CVP_POLAR = "NQST"
# Wilkinson–Harrison çözünürlük modelinde dönüş oluşturan kalıntılar.
TURN_FORMING_RESIDUES = "NGPS"

# TOP-IDP düzensizlik eğilimi ölçeği (Campen ve ark. 2008); pozitif değer
# düzensizliği, negatif değer düzenli yapıyı destekler.
TOP_IDP_SCALE = {
    "A": 0.06, "C": 0.02, "D": 0.192, "E": 0.736, "F": -0.697, "G": 0.166, "H": 0.303,
    "I": -0.486, "K": 0.586, "L": -0.326, "M": -0.397, "N": 0.007, "P": 0.987, "Q": 0.318,
    "R": 0.180, "S": 0.341, "T": 0.059, "V": -0.121, "W": -0.884, "Y": -0.510,
}


# tAI wobble cezaları s (dos Reis ve ark. 2004, S. cerevisiae için optimize):
# "antikodon 34. baz : kodon 3. baz". I inozin (A34 deaminasyonu), L lizidin
# (prokaryotlarda tRNA-Ile CAU'nun AUA'yı okuyan biçimi).
TAI_S_VALUES = {
    "I:U": 0.0, "G:C": 0.0, "U:A": 0.0, "C:G": 0.0,
    "G:U": 0.41, "I:C": 0.28, "I:A": 0.9999, "U:G": 0.68, "L:A": 0.89,
}
# Kodonun 3. bazını okuyabilen antikodonlar: kodon 3. baz -> ((antikodonun
# 34. bazı genomda yazıldığı gibi, TAI_S_VALUES anahtarı), ...). Genomdaki A34
# olgun tRNA'da I'ya dönüşür.
TAI_WOBBLE_PAIRS = {
    "U": (("A", "I:U"), ("G", "G:U")),
    "C": (("G", "G:C"), ("A", "I:C")),
    "A": (("U", "U:A"), ("A", "I:A")),
    "G": (("C", "C:G"), ("U", "U:G")),
}
