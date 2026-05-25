# Audit Komputasi Perselisihan Efisiensi Proteom 2024: Efek Warburg

Repositori ini berisi kode sumber untuk audit komputasi terhadap dua makalah tahun 2024 yang menyimpulkan hal berlawanan tentang efisiensi proteom jalur produksi ATP:

- **Shen dkk.** (*Nat. Chem. Biol.* 2024) — respirasi lebih efisien per massa protein
- **Kukurugya dkk.** (*PNAS* 2024) — glikolisis menghasilkan ATP lebih cepat per miligram enzim

Pipeline ini mengimplementasikan ulang kedua model dalam satu basis kode, menjalankan tiga eksperimen diagnostik, dan mengidentifikasi pengukuran penuntas yang akan menyelesaikan perdebatan.

## Struktur Proyek

```
├── main.py                  # Entry point — jalankan seluruh pipeline
├── src/
│   ├── params.py            # M2: Parameter provenance dari literatur
│   ├── models.py            # M1: Solver (LP Shen + analitik Kukurugya)
│   ├── audit.py             # M3: Eksperimen E1, E2, E2c, E2d
│   ├── identifiability.py   # M4: Analisis sensitivitas Sobol + Morris
│   ├── viz.py               # M5: Visualisasi (6 figur)
│   └── eccore_validation.py # M6: Validasi ecCore FBA (opsional)
├── tests/
│   ├── test_models.py       # Unit test modul M1–M4
│   └── test_system.py       # System test end-to-end
├── results/                 # Output CSV dari pipeline
├── figures/                 # Output figur (PNG + SVG)
├── doc/                     # Laporan LaTeX (format IEEE)
└── requirements.txt
```

## Instalasi

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Penggunaan

Jalankan seluruh pipeline (M2 → M1 → M3 → M4 → M5):

```bash
python main.py
```

Untuk menyertakan validasi ecCore FBA (M6, membutuhkan COBRApy):

```bash
python main.py --stretch
```

Hasil akan ditulis ke `results/` (CSV) dan `figures/` (PNG + SVG).

## Pengujian

```bash
pytest
```

## Pipeline

| Modul | Deskripsi |
|-------|-----------|
| **M2** | Ekstraksi dan dokumentasi provenance parameter dari SI ketiga paper |
| **M1** | Dua solver dengan antarmuka bersama: LP 2-sektor (Shen) dan analitik 5-parameter (Kukurugya) |
| **M3** | Eksperimen audit: overlay sumbu-bersama (E1), kapasitas vs. realisasi (E2), bootstrap CI (E2c), dekomposisi atribusi (E2d) |
| **M4** | Analisis sensitivitas global: Sobol (S₁, Sₜ) dan Morris (μ*, σ) |
| **M5** | Visualisasi seluruh hasil ke 6 figur publikasi |
| **M6** | Validasi independen dengan model FBA berbatas-enzim ecCore (opsional) |

## Temuan Utama

- Perselisihan bersumber pada **besaran yang diukur**: kapasitas V·γ (Kukurugya) vs. efisiensi terealisasi u·V·γ (Shen)
- Verdict membalik pada u_G ≈ ρ untuk ketiga organisme (90% CI seluruhnya < 1.0)
- Analisis Sobol mengidentifikasi **u_G** sebagai pengukuran penuntas (Sₜ ≈ 0.76–0.78)
- Mekanisme penyilangan komplementer terhadap pendamaian Wang 2025

## Kesimpulan

Shen dan Kukurugya tidak berselisih tentang biologi. Keduanya benar, hanya mengukur hal yang berbeda. Kukurugya mengukur **kapasitas maksimal** enzim (V·γ, seolah semua enzim bekerja penuh), sementara Shen mengukur **efisiensi terealisasi** (u·V·γ, memperhitungkan bahwa sebagian enzim glikolitik menganggur di dalam sel). Ketika enzim glikolitik tidak bekerja 100%, keunggulan glikolisis menyusut, dan pada titik tertentu (u_G ≈ ρ), respirasi menjadi lebih efisien. Itulah yang diamati Shen.

Satu angka yang akan menyelesaikan perdebatan ini selamanya: **u_G** — berapa persen kapasitas enzim glikolitik yang benar-benar dipakai sel secara in-vivo. Ukur itu, dan pertanyaan "jalur mana yang lebih efisien?" langsung terjawab.

## Penulis

- Naufarrel Zhafif Abhista (13523149)
- Frederiko Eldad Mugiyono (13523147)
- Hasri Fayadh Muqaffa (13523156)
- I Made Wiweka Putera (13523160)

Program Studi Teknik Informatika, Institut Teknologi Bandung
