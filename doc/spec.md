# IF3211 — Cetak Biru Proyek
## Topik 6 (Respirasi Sel & Produksi ATP) · Soal 7: Metabolisme Berbatas-Proteom & Efek Warburg

> **Kalimat pemikat:** Literatur 2024 memuat *kontradiksi langsung yang belum terselesaikan* tentang apakah respirasi atau glikolisis lebih "efisien-proteom". Kami **tidak** membangun model baru untuk berpihak. Kami membangun satu *harness* (kerangka uji) sumber-terbuka, bebas-AI, yang **mengimplementasikan ulang kedua model bersaing 2024 secara berdampingan**, memberi mereka input identik, lalu **mendiagnosis dari mana perbedaannya berasal** — dan menentukan **eksperimen tunggal mana yang akan menyelesaikan perdebatan itu.**

---

## 0. TL;DR (baca ini dulu)

- **Apa:** Sebuah *audit komputasi terkendali* terhadap perselisihan efisiensi-proteom 2024. Kami mengimplementasikan ulang **model 2-sektor Shen et al.** dan **model 5-parameter Kukurugya et al.** dalam satu basis kode dengan antarmuka input bersama, menjalankan keduanya pada input identik, lalu menambahkan dua eksperimen diagnostik: (1) *batas akuntansi* protein, (2) *analisis identifiabilitas* (sensitivitas global).
- **Mengapa ini baru (nilai jual):** Wang (eLife 2025) sudah "mendamaikan" kedua makalah dengan *menambah model ketiga* (alokasi protein + heterogenitas sel). Belum ada yang melakukan **perbandingan kepala-ke-kepala terhadap dua model asli** untuk melokalisasi *sumber* perbedaan (definisi vs. parameter vs. struktur), dan belum ada yang menanyakan **pengukuran mana yang akan membalikkan kesimpulan** (kerangka "nilai-informasi"). Itu kontribusi meta-ilmiah yang jujur, dapat dipertahankan, dan layak dikerjakan mahasiswa S1 dalam 5 hari.
- **Bukan AI/ML:** Murni optimasi matematis (LP via `scipy`/`PuLP`-CBC) + analisis sensitivitas Monte Carlo klasik (Sobol/Morris via SALib). Sengaja demikian — transparansi analitis adalah inti argumen kami.
- **Tumpukan teknologi:** Python ≥3.10, NumPy, SciPy, PuLP (solver CBC — sumber terbuka, **tidak perlu lisensi Gurobi**), SALib, Matplotlib, pandas. Opsional (peregangan): COBRApy untuk model inti berbatas-enzim.
- **Luaran:** (1) Basis kode GitHub, sekali-jalan dan dapat direproduksi; (2) **paper format IEEE ≤6 halaman** (ditulis dalam Markdown berstruktur IEEE, dua kolom saat dikompilasi); (3) video ≤10 menit dengan dek slide. Spesifikasi ketiganya ada di §9–§11.
- **Disiplin lingkup:** Inti (wajib-kirim) = Lapisan 1–4. Peregangan (hanya jika lebih cepat dari jadwal) = Lapisan 5. Urutan pemangkasan didefinisikan di §12.

---

## 1. Pertanyaan Penelitian & Hipotesis

**Ketegangan latar.** Dua makalah 2024 yang sama-sama tinjauan-sejawat menyimpulkan hal berlawanan; sebuah makalah 2025 mencoba mendamaikan dengan model *baru*:

| Sumber | Klaim | Arah (`ρ = ε_R/ε_G`) |
|---|---|---|
| Shen et al. 2024 (*Nat. Chem. Biol.* 20:1123) | Respirasi "beberapa kali lipat lebih efisien-proteom" daripada glikolisis | `ρ > 1` |
| Kukurugya, Rosset & Titov 2024 (*PNAS* 121(46):e2409509121) | Glikolisis menghasilkan ATP **0,54× / 2,1× / 3,1×** lebih cepat per mg enzim daripada respirasi pada *E. coli* / ragi / sel mamalia; model 5-parameter | E. coli `ρ≈1,85`, ragi `ρ≈0,48`, mamalia `ρ≈0,32` |
| Wang 2025 (*eLife* 13:RP94586) | **Model alokasi-protein + heterogenitas sel baru** yang mendamaikan keduanya: efisiensi respirasi & fermentasi berpotongan menurut kualitas nutrisi | mendamaikan via penyilangan `ρ` |

> **⚠ Peringatan integritas (wajib dibaca).** Angka `ρ` di tabel di atas adalah **provisional**, diturunkan dengan membalik kalimat ringkas — **bukan dikutip dari tabel parameter asli.** Sebelum dipakai di paper, **wajib** diverifikasi dari tabel/SI masing-masing makalah (lihat gerbang Hari-0 di §8). Khususnya: Kukurugya menyatakan **E. coli adalah kasus khusus** (memakai jalur respiro-fermentatif Pta-AckA, *tidak pernah* beralih ke fermentasi murni walau O₂ ada) — jadi jangan klaim "E. coli `ρ≈1,85`" tanpa verifikasi; tandanya bisa berbeda.
>
> **Kehalusan yang sebenarnya kami eksploitasi (sumber perselisihan yang lebih mungkin):** Kukurugya mengukur **kapasitas maksimal** efisiensi, yaitu `V·γ` dengan `V` = aktivitas spesifik maksimal (µmol substrat·menit⁻¹·mg protein⁻¹, alias Vmax). Shen mengukur **efisiensi ter-realisasi in-vivo**, yaitu fluks ATP aktual ÷ massa enzim terukur. Keduanya **menyebutnya** "proteome efficiency" (Kukurugya bahkan mengutip Shen sebagai metrik yang sama), padahal **kapasitas ≠ realisasi**. Shen sendiri menemukan protein glikolitik diekspresikan **berlebih/menganggur** ("constitutive glycolysis"), sehingga efisiensi terealisasinya rendah meski kapasitas per proteinnya tinggi. Inilah dugaan akar perselisihan — bukan sekadar "enzim mana yang dihitung". Audit kami menguji dugaan ini secara langsung (lihat E2 di §3.3).

**Mengapa kami TIDAK sekadar "mendamaikan".** Pendamaian sudah diterbitkan (Wang 2025). Mengulanginya bukan kontribusi. Sebaliknya, kami menempatkan diri sebagai **auditor**: menjalankan kembali kedua model asli pada pijakan yang sama untuk mendiagnosis akar perselisihan dan eksperimen yang dapat menuntaskannya.

**Pertanyaan penelitian utama (RQ1).**
> Bila model Shen dan model 5-parameter Kukurugya dijalankan pada **input identik** dalam satu kerangka, apakah keduanya benar-benar bertentangan? Dan dari empat kemungkinan sumber perbedaan ini — **(i) besaran yang diukur** (kapasitas maksimal `V·γ` vs. fluks ter-realisasi `J/m_enzim`), **(ii) definisi** "protein jalur" (enzim mana yang dihitung), **(iii) nilai parameter** terukur, **(iv) struktur model** — manakah yang menjadi sumber dominan?

**Pertanyaan penelitian kedua (RQ2 — kerangka nilai-informasi).**
> Parameter terukur tunggal manakah yang, bila dipastikan, **membalikkan kesimpulan** "respirasi vs. glikolisis lebih efisien"? Dengan kata lain: **eksperimen mana yang akan menyelesaikan perdebatan ini?**

**Hipotesis (H).**
> Perselisihan Shen/Kukurugya **bukan** artefak struktur model, melainkan terutama karena **kedua makalah mengukur besaran yang berbeda**: Kukurugya memakai kapasitas maksimal `V·γ` sedangkan Shen memakai fluks ter-realisasi `J/m_enzim`. Kami berhipotesis bahwa **eksperimen E2 (kapasitas vs. realisasi)** akan menunjukkan verdict membalik ketika fraksi enzim menganggur (idle) dinaikkan dari 0 ke nilai realistis, dan bahwa **analisis identifiabilitas E3** akan menempatkan parameter rasio kapasitas (`V·γ` glikolisis vs. respirasi) — atau fraksi idle — sebagai pengendali utama verdict. Sumbu "definisi enzim" diuji sebagai hipotesis sekunder. Hasilnya menamai **pengukuran yang akan menuntaskan perdebatan**.

Pasangan RQ/H ini adalah tulang punggung laporan maupun dek.

---

## 2. Konsep Biologi (yang harus Anda tunjukkan paham)

Singkat saja — ini 1 halaman laporan, bukan inti proyek.

- **Respirasi vs. fermentasi.** Glukosa → 2 piruvat (glikolisis, sitosol, ~2 ATP neto, tak bergantung O₂). Piruvat → CO₂ via siklus TCA + fosforilasi oksidatif (ETC + ATP sintase, mitokondria, ~30–32 ATP, bergantung O₂). Respirasi memiliki **hasil ATP per glukosa** jauh lebih tinggi tetapi memakai aparatus enzim yang besar dan lambat.
- **Efek Warburg.** Banyak sel tumbuh-cepat (tumor; ragi — "efek Crabtree") memfermentasi glukosa menjadi laktat/etanol *walau O₂ tersedia*. Berlawanan intuisi, sebab "membuang" ~15× hasil ATP.
- **Penjelasan efisiensi-proteom.** Hipotesis klasik (Basan 2015; Pfeiffer 2001): sel punya **anggaran protein terbatas**; enzim glikolitik diduga memberi lebih banyak ATP per gram enzim per jam daripada rantai pernapasan yang besar, sehingga di bawah anggaran ketat, glikolisis aerobik memaksimalkan laju ATP.
- **Empat besaran yang WAJIB dibedakan (jangan dicampur — ini sumber kebingungan literatur).** (1) **Yield** `γ` = mol ATP per mol glukosa (respirasi ≫ glikolisis). (2) **Aktivitas spesifik / kapasitas** `V` = µmol substrat·menit⁻¹·mg protein⁻¹ pada saturasi (Vmax). (3) **Laju ATP per protein (kapasitas)** = `V·γ` — *inilah* yang Kukurugya sebut "lebih cepat per gram protein". (4) **Efisiensi ter-realisasi** = fluks ATP **aktual in-vivo** ÷ massa enzim **terukur** — inilah yang Shen ukur. Besaran (3) dan (4) **berbeda**: enzim bisa berkapasitas tinggi per protein namun beroperasi jauh di bawah kapasitas (idle), sehingga realisasinya rendah. Membedakan keempatnya adalah inti pemahaman biologi proyek ini.
- **Pertarungan 2024** (lihat §1). Shen: respirasi lebih efisien **secara terealisasi** → protein glikolitik diekspresikan berlebih sebagai cadangan ("*proteome hedging*" / robust terhadap hipoksia). Kukurugya: glikolisis lebih cepat **secara kapasitas** `V·γ` → beralih saat glukosa berlimpah & proteom jenuh. Wang (2025): efisiensi respirasi & fermentasi **berpotongan** menurut kualitas nutrisi + heterogenitas sel.

**Yang Anda modelkan:** sel mengalokasikan anggaran protein tetap di antara dua sektor penghasil ATP untuk memaksimalkan **laju** ATP, dengan **ketersediaan glukosa + anggaran proteom** sebagai kendala pengikat (**bukan** O₂ — Warburg terjadi justru saat O₂ ada). Audit kami membandingkan bagaimana dua makalah memetakan biologi yang sama ke besaran efisiensi yang berbeda.

---

## 3. Model & Eksperimen Audit (spesifikasi matematis)

Ini seluruh inti intelektual. Sengaja kecil dan transparan.

### 3.1 Model A — alokasi proteom berbasis-glukosa (gaya Shen/Basan)
Variabel: fraksi massa proteom `φ_G` (glikolisis), `φ_R` (respirasi), `φ_Q` tetap (housekeeping). Parameter selaras notasi Kukurugya agar dapat dibandingkan: `V_i` (aktivitas spesifik/kapasitas), `γ_i` (yield ATP/glukosa), `Φ = φ_G+φ_R` (fraksi proteom maksimal untuk enzim penghasil-ATP).
```
Laju serapan glukosa per sektor:  J_glc,i ≤ V_i · φ_i        (i ∈ {G,R})
Laju ATP per sektor:              J_ATP,i = γ_i · J_glc,i
(C1) Anggaran proteom:            φ_G + φ_R ≤ Φ
(C2) Ketersediaan glukosa:        J_glc,G + J_glc,R ≤ g_avail   ← KNOB UTAMA (bukan O₂)
(C3) Non-negatif:                 φ_G, φ_R, J_glc,i ≥ 0
Tujuan:  maksimumkan J_ATP = γ_G·J_glc,G + γ_R·J_glc,R
```
Ini **LP**; optimum di titik-sudut → peralihan respirasi↔glikolisis yang bersih saat `g_avail` naik. **Catatan koreksi:** versi spec sebelumnya memakai kapasitas O₂ sebagai pengikat — itu keliru secara biologi (glikolisis aerobik = O₂ ada). Knob yang benar adalah **ketersediaan glukosa**.

### 3.2 Model B — model 5-parameter Kukurugya (analitik, tanpa parameter bebas)
Lima parameter: **γ_resp, γ_glyc** (yield ATP/glukosa), **V_resp, V_glyc** (aktivitas spesifik = laju substrat per mg protein-jalur), **Φ** (fraksi proteom maksimal untuk enzim penghasil-ATP). Model memaksimalkan laju ATP terhadap alokasi proteom pada tiap tingkat ketersediaan glukosa dan memiliki **solusi analitik unik tanpa parameter yang dapat disetel**. Aturan inti yang harus direproduksi: pada glukosa rendah, respirasi menang (yield tinggi `γ_resp>γ_glyc`); pada glukosa tinggi & proteom jenuh, jalur dengan **laju ATP per protein lebih tinggi** (`V·γ` terbesar) menang. Untuk *E. coli* model **tidak** memprediksi peralihan ke fermentasi murni (kasus Pta-AckA). **Lima parameter wajib diangkat dari tabel/SI PNAS — bukan ditebak** (gerbang Hari-0, §8). Implementasikan **apa adanya**, tanpa "memperbaiki".

### 3.3 Eksperimen audit (kontribusi baru)

**E1 — Harness sumbu-bersama (replikasi).**
Beri Model A dan Model B *vektor input identik* dan plot keluaran keduanya pada **sumbu observabel yang sama**: sumbu-x = **laju serapan glukosa** `g_avail`; sumbu-y = **fraksi fluks fermentatif** (dan/atau laju ATP). Pertanyaan: di mana ambang peralihan keduanya berbeda, dan apakah perbedaan itu kuantitatif (lokasi ambang) atau kualitatif (arah verdict)? Observabel bersama ini **wajib didefinisikan eksplisit di kode**; tanpa itu "overlay" tak bermakna.

**E2 — Kapasitas vs. realisasi (uji utama "besaran yang diukur").**
Inilah diagnosis inti yang dikoreksi. Definisikan **fraksi pemanfaatan kapasitas** `u_i ∈ (0,1]` per sektor (1 = beroperasi di Vmax; <1 = ada enzim idle/cadangan). Hitung dua efisiensi pada parameter yang sama: **kapasitas** `(V·γ)_i` (gaya Kukurugya, `u=1`) dan **terealisasi** `u_i·(V·γ)_i` (gaya Shen). Sapu `u_G` (mencerminkan ekspresi-berlebih glikolitik yang dilaporkan Shen) dan **temukan nilai `u_G` di mana verdict `(V·γ)_R` vs. `u_G·(V·γ)_G` membalik.** Bila verdict membalik pada `u_G` realistis, perselisihan terjelaskan oleh kapasitas-vs-realisasi — bukan kesalahan salah satu pihak. *Sub-eksperimen sekunder:* variasikan atribusi enzim membran/pemeliharaan (sumbu "definisi") dan laporkan apakah ini menambah pembalikan di luar efek `u_G`.

**E3 — Identifiabilitas / nilai-informasi (uji "parameter"), well-posed.**
Jalankan sensitivitas global (Sobol + Morris via **SALib** — Monte Carlo, *bukan* ML) dengan **output KONTINU**, bukan biner: gunakan **margin efisiensi** `m = u_R·(V·γ)_R − u_G·(V·γ)_G` (atau ambang glukosa kritis `g*` tempat peralihan terjadi). Distribusi input parameter diambil dari **rentang ketidakpastian di SI tiap makalah** (bukan rentang sembarang). Peringkatkan parameter menurut indeks Sobol total-order → parameter teratas = **pengukuran yang, bila dipastikan, paling mengurangi ketidakpastian verdict** = eksperimen penuntas. (Catatan posisi & landasan: kerangka sensitivitas untuk model alokasi-protein, *sEnz* (Bioinformatics 2024), memakai koefisien sensitivitas kapasitas/enzim via shadow price; kami memakai Sobol/Morris yang lebih sederhana dan **menerapkannya spesifik pada perselisihan Shen–Kukurugya**, yang belum dilakukan.)

**Gerbang validasi (WAJIB sebelum audit apa pun).** Tiap model hasil re-implementasi **hanya boleh masuk audit** setelah mereproduksi ≥1 hasil terbit dalam toleransi tertera: Model B → bentuk kurva peralihan & laju glikolisis/respirasi vs. glukosa (Kukurugya Fig. 1F–H) dalam ±20%; Model A → perilaku batas (glukosa→0 ⇒ respirasi murni; proteom jenuh & glukosa tinggi ⇒ onset fermentasi). Tanpa lolos gerbang ini, perbandingan tidak sah secara akademik.

**E3b — Sobol pada keluaran LP nonlinier (uji non-tautologis).**
E3 asli menjalankan Sobol pada margin `m = u_R·(V·γ)_R − u_G·(V·γ)_G` — rumus yang **linier** di `u_G`, sehingga Sobol secara trivial menempatkan `u_G` di peringkat pertama (tautologis). E3b menjalankan Sobol pada **keluaran LP Model A** (`frac_glyc` = fraksi fermentatif optimal) yang **nonlinier** — memiliki transisi titik-sudut dan pergantian kendala aktif. Bila `u_G` masih mendominasi pada keluaran nonlinier, temuan E3 bersifat robust dan bukan artefak rumus linier. Glukosa di-set pada `g* × 1.5` (sedikit di atas ambang transisi) agar LP beroperasi di rezim di mana keputusan fenotipe benar-benar terjadi.

**E4 — Diagram fase 2D (peta rezim metabolik).**
Hitung fenotipe LP pada grid 2D `(u_G, g_avail)` per organisme. Visualisasikan sebagai heatmap kontinu (fraksi fermentatif) dengan kontur batas rezim (respirasi / campuran / glikolisis). Diagram ini **belum ada** di literatur Shen, Kukurugya, maupun Wang dan merupakan artefak visual orisinal yang menunjukkan bagaimana kedua sumbu (pemanfaatan kapasitas dan ketersediaan glukosa) berinteraksi secara simultan untuk menentukan fenotipe metabolik.

**Punchline yang dapat difalsifikasi (hasil utama yang diharapkan):**
Perselisihan terlokalisasi terutama pada **besaran yang diukur (kapasitas vs. realisasi)**, bukan struktur model; verdict membalik pada fraksi-idle glikolitik `u_G` yang realistis (konsisten dengan temuan ekspresi-berlebih Shen); dan analisis identifiabilitas menamai pengukuran tunggal yang menuntaskan perdebatan. Posisi terhadap Wang 2025: ia *menambah model* (penyilangan efisiensi + heterogenitas) untuk mendamaikan; kami *mengaudit model asli* untuk mendiagnosis bahwa keduanya mungkin tak mengukur besaran yang sama. Diagnosis berbeda, komplementer — dan dapat dibandingkan langsung dengan mekanisme Wang.

---

## 4. Arsitektur Sistem

```
                ┌──────────────────────────────────────────────┐
                │                main.py / run.ipynb           │
                │   (satu titik masuk → regenerasi semua gbr)  │
                └───────────────┬──────────────────────────────┘
                                │
     ┌────────────┬─────────────┴───────────┬───────────────┬─────────────┐
     ▼            ▼                         ▼               ▼             ▼
┌──────────┐  ┌──────────────┐      ┌─────────────────┐ ┌──────────┐  ┌────────┐
│ M2 params│  │ M1 models    │      │ M3 audit        │ │ M4 ident.│  │ M5 viz │
│ ε dari 3 │─▶│ A: 2-sektor  │─────▶│ E1 input-bersama│ │ E3 Sobol/│─▶│ gambar │
│ makalah  │  │ B: 5-param   │      │ E2 batas-       │ │ Morris   │  │ utk    │
│+provenans│  │ antarmuka    │      │    akuntansi    │ │ (SALib)  │  │ lap/dek│
└──────────┘  │ bersama      │      └─────────────────┘ └──────────┘  └────────┘
              └──────┬───────┘
                     ▼ (peregangan/validasi saja)
              ┌──────────────────┐
              │ M6 ecCore (COBRA)│  inti berbatas-enzim, opsional
              └──────────────────┘
```

**Aliran data:** `M2` memasok parameter ketiga makalah → `M1` solver murni (dua model, satu antarmuka) → `M3` menjalankan E1 & E2 → `M4` menjalankan E3 (sensitivitas) → `M5` merender. Semua deterministik kecuali sampling Monte Carlo M4 (yang diberi seed tetap untuk reproduksibilitas).

---

## 5. Spesifikasi Per-Modul

Tiap modul = satu berkas Python, satu pemilik, satu uji-terima jelas. "Tujuan unit" = pemeriksaan yang dapat difalsifikasi bahwa unit selesai.

### M1 — `src/models.py` (dua solver, satu antarmuka)
- **Tugas unit:** Implementasikan `solve_shen(params)` (LP 2-sektor) dan `solve_kukurugya(params)` (model 5-parameter), keduanya menerima `dataclass ModelParams` yang sama dan mengembalikan `ModelResult` yang sama (`phenotype`, `eps_ratio`, `J_ATP`, `frac_glyc`).
- **Peran sistem:** Satu-satunya tempat matematika model tinggal; dipanggil di mana-mana.
- **Implementasi:** `scipy.optimize.linprog` (HiGHS) untuk Model A; aritmetika tertutup/akar untuk Model B sesuai PNAS.
- **Tujuan unit:** Pada O₂ tak terbatas & `ρ>1`, Model A → respirasi murni; O₂→0 → glikolisis murni. Model B mereproduksi rasio 0,54/2,1/3,1 Kukurugya untuk E. coli/ragi/mamalia hingga toleransi yang dilaporkan.
- **Pemilik:** Anggota A. **Estimasi:** 1 hari (B perlu pembacaan PNAS cermat).

### M2 — `src/params.py` (provenans parameter)
- **Tugas unit:** Kodekan `ε_G`, `ε_R`, hasil, kapasitas untuk **ketiga makalah** × tiga organisme, tiap angka ditandai sitasi + nomor gambar/tabel/SI.
- **Tujuan unit:** `python -m src.params --table` mencetak tabel provenans bersih; setiap angka tertelusur. Tidak ada "konstanta sihir".
- **Pemilik:** Anggota B. **Estimasi:** 0,5 hari (terutama membaca SI dengan teliti).

### M3 — `src/audit.py` (E1 + E2 — inti baru, bagian 1)
- **Tugas unit:** E1: jalankan M1.A & M1.B pada grid input identik, kumpulkan keluaran ke satu dataframe pada sumbu bersama. E2: variasikan atribusi enzim (batas akuntansi), hitung `ε_R/ε_G(atribusi)`, temukan titik balik verdict.
- **Tujuan unit:** Menghasilkan `audit_overlay.parquet` dan `accounting_flip.csv`; titik-balik verdict terdeteksi & dilaporkan numerik.
- **Pemilik:** Anggota C. **Estimasi:** 1 hari.

### M4 — `src/identifiability.py` (E3 — inti baru, bagian 2)
- **Tugas unit:** Definisikan ruang parameter + rentang; jalankan Sobol & Morris (SALib) dengan keluaran = indikator verdict; peringkatkan parameter; identifikasi "pengukuran penuntas".
- **Tujuan unit:** Menghasilkan `sobol_indices.csv` + peringkat; parameter teratas konsisten antara Sobol & Morris (pemeriksaan silang metode).
- **Pemilik:** Anggota C atau D. **Estimasi:** 1 hari (alur SALib standar).

### M5 — `src/viz.py` (gambar)
- **Tugas unit:** Render: (F1) overlay E1 dua model pada sumbu bersama; (F2) kurva balik batas-akuntansi E2; (F3) tornado/diagram indeks Sobol E3 ("pengukuran penuntas"); (F4 opsional) peralihan ecCore.
- **Tujuan unit:** `python -m src.viz` meregenerasi semua gambar laporan dari keluaran tercache dalam <30 dtk.
- **Pemilik:** Anggota D. **Estimasi:** 1 hari (tumpang-tindih M3/M4).

### M6 — `src/eccore_validation.py` (PEREGANGAN — validasi mekanistik)
- **Tugas unit:** Muat model inti terbit (`e_coli_core` via COBRApy), tambah satu kendala massa-enzim `Σ mᵢ·vᵢ ≤ Φ`, sapu `Φ`/glukosa, konfirmasi peralihan respirasi→fermentasi yang diprediksi model mainan.
- **Tujuan unit:** ecFBA menunjukkan onset limpahan asetat saat `Φ` mengetat, secara kualitatif cocok arah audit.
- **Pemilik:** Anggota A/C, **hanya bila M1–M5 terkirim pada Hari 3.** **Estimasi:** 1 hari. **Dipangkas pertama bila terlambat.**

### Orkestrasi — `main.py`
- **Tugas unit:** `python main.py` menjalankan M2→M1→M3→M4→M5 (dan M6 bila `--stretch`), menulis semua artefak, mencetak temuan numerik utama (sumber perbedaan; parameter teratas).
- **Tujuan unit:** Klon segar + `pip install -r requirements.txt` + `python main.py` mereproduksi setiap gambar tanpa langkah manual.

---

## 6. Tumpukan Teknologi (eksplisit non-AI/ML)

| Lapisan | Pilihan | Alasan |
|---|---|---|
| Bahasa | Python ≥ 3.10 | wajib oleh spesifikasi |
| Solver LP | `scipy.optimize.linprog` (HiGHS); cadangan `PuLP` + CBC | sumber terbuka, tanpa lisensi, tanpa AI |
| Sensitivitas | **SALib** (Sobol, Morris) | Monte Carlo klasik; metode statistik, bukan ML |
| Numerik | NumPy | grid, sapuan tervektorisasi |
| Data | pandas + parquet/CSV | provenans & reproduksibilitas |
| Plot | Matplotlib (+ opsional seaborn) | gambar kualitas publikasi |
| GEM peregangan | COBRApy + `e_coli_core` | validasi mekanistik |
| Reproduksi | `requirements.txt` versi terkunci, `main.py`, seed tetap | sekali-jalan |
| VCS | Git + repo GitHub yang dapat diakses dosen+asisten | wajib oleh spesifikasi |

**Sengaja dikecualikan:** semua jaring saraf, regresi-sebagai-model, surrogate ML, AutoML. Justifikasi (nyatakan dalam laporan): nilai proyek ini adalah *transparansi analitis* — LP bertitik-sudut yang bermakna biologis langsung, plus sensitivitas Sobol yang dapat diinterpretasi. ML justru mengaburkan interpretabilitas yang membuat audit ini meyakinkan. Tidak perlu Gurobi/CPLEX: masalahnya kecil; CBC/GLPK/HiGHS menanganinya seketika.

---

## 7. Struktur Repositori

```
warburg-audit/
├── README.md                 # apa/mengapa/cara-jalan, galeri gambar, tautan video
├── requirements.txt          # versi terkunci
├── LICENSE                   # MIT
├── main.py                   # pipeline sekali-jalan
├── src/
│   ├── models.py             # M1 (Model A & B, antarmuka bersama)
│   ├── params.py             # M2
│   ├── audit.py              # M3 (E1 + E2)
│   ├── identifiability.py    # M4 (E3, SALib)
│   ├── eccore_validation.py  # M6 (peregangan)
│   └── viz.py                # M5
├── data/
│   └── provenance.csv        # sumber parameter (otomatis)
├── results/                  # keluaran audit & sensitivitas
├── figures/                  # F1–F7 png+svg
├── tests/
│   └── test_models.py        # uji solusi-sudut + reproduksi rasio Kukurugya
└── report/                   # PDF 6-halaman + dek slide
```

---

## 8. Linimasa — Sprint 5 Hari (24 Mei → 29 Mei 2026)

Realistis untuk 3–4 orang bekerja paralel. **Catatan kritis: hari ini 24 Mei, batas akhir 29 Mei — ini sprint, bukan proyek 14 minggu.**

| Hari | Tanggal | Tujuan | Pemilik | Selesai bila… |
|---|---|---|---|---|
| **0** | Min 24 Mei (pagi) | **GERBANG PARAMETER & KONSEP.** Angkat 5 parameter Kukurugya (γ_resp, γ_glyc, V_resp, V_glyc, Φ) dari tabel/SI PNAS + verifikasi arah `ρ` per organisme; pastikan seluruh tim paham beda 4 besaran (γ, V, V·γ, efisiensi terealisasi) | B (SI), semua (konsep) | Tabel 5 parameter terisi dari sumber asli; tak ada angka tebakan |
| **1** | Min 24 Mei | Kerangka repo; Model A (LP berbasis-glukosa) jalan; mulai Model B (analitik); M2 provenans | A (Model A), B (Model B+M2), C+D (repo+pustaka) | Model A lulus uji batas (glukosa→0 ⇒ respirasi murni) |
| **2** | Sen 25 Mei | Model B jalan & **lolos gerbang validasi** (reproduksi Kukurugya Fig.1F–H ±20%); E1 overlay sumbu-bersama; mulai paper §I–§II | A/B (Model B+validasi), C (E1), B (paper) | F1 overlay dua-model ada; gerbang validasi lolos |
| **3** | Sel 26 Mei | E2 kapasitas-vs-realisasi (F2, sapu fraksi idle `u_G`); mulai E3 Sobol (output kontinu); **putuskan lampu hijau/merah peregangan** | C (E2), D (E3), A (uji) | F2 ada (titik-balik `u_G` terdeteksi); E3 berjalan |
| **4** | Rab 27 Mei | Selesaikan E3 (F3 "pengukuran penuntas"); (peregangan M6 bila hijau); tulis **paper §III–§VI**; bersihkan kode+uji+README | C/D (E3+viz), B (paper), A (M6/uji) | Draf paper lengkap; `python main.py` reproduksi semua |
| **5** | Kam 28 Mei | Bangun dek; rekam video ≤10 mnt (semua anggota tampil); poles **paper final** (cek notasi, referensi, anggaran halaman); dorong repo, verifikasi akses incognito | semua | Video terunggah, tautan di paper; repo dapat diakses dosen+asisten |
| **kirim** | **Jum 29 Mei** | Kirim **paper** + kode + berkas + tautan video | — | sesuai kanal dosen |

**Standup harian (15 mnt):** apa yang terkirim, apa yang macet, apakah garis-pangkas (§12) terpicu.

---

## 9. Luaran 1 — Spesifikasi Basis Kode

- **Aksesibilitas:** GitHub publik (atau dibagikan ke email penilai). Uji tautan di mode incognito sebelum pengumpulan — ini butir penilaian eksplisit.
- **Standar reproduksibilitas:** klon segar → `pip install -r requirements.txt` → `python main.py` meregenerasi tiap gambar. Tanpa keadaan notebook tersembunyi, tanpa urutan sel manual. Seed tetap untuk M4.
- **README wajib memuat:** kalimat-pemikat, RQ, cara jalan, galeri gambar (sematkan F1–F3), catatan provenans parameter, tautan video, tabel tim + kontribusi, **tautan ke paper (PDF)**.
- **Uji:** minimal `tests/test_models.py` (solusi-sudut Model A; reproduksi rasio Kukurugya Model B). Jalankan `pytest`. Asuransi murah terhadap "tidak jalan" saat dinilai.
- **Higiene:** type hints + docstring pada M1–M4 (penilai membacanya), tanpa kode mati terkomentari, `results/` dan `figures/` di-commit agar penilai tak perlu menjalankan apa pun untuk melihat keluaran.

---

## 10. Luaran 2 — Spesifikasi Paper (Format IEEE, ≤ 6 halaman)

**Format:** IEEE two-column conference style. Ditulis dalam Markdown dengan heading `##` sebagai penanda section IEEE (I, II, III, …). Untuk tampilan final dua-kolom, render via **Pandoc + template IEEEtran** (`pandoc paper.md -o paper.pdf --template ieee`) atau paste ke **Overleaf** dengan template IEEEtran. Kalau tidak ada yang bisa LaTeX, Word template IEEE (tersedia di ieee.org) juga diterima — yang penting strukturnya.

> **Konfirmasi ke dosen:** pastikan "format bebas, maksimal 6 halaman" berlaku juga untuk format IEEE. IEEE two-column 6 halaman itu sangat lega (~7.000 kata) — kemungkinan besar aman, tapi tanyakan lebih dulu.

---

### Struktur paper & anggaran halaman

| No. | Section IEEE | Setara spesifikasi tugas | Hlm (est.) | Isi |
|---|---|---|---|---|
| — | **Judul, Penulis, Afiliasi** | — | — | Judul dalam bahasa Indonesia + Inggris; nama + NIM semua anggota; institusi; email ketua kelompok. |
| — | **Abstrak** | — | 0,15 | Dua paragraf, ≤250 kata. Para. 1: konteks (perdebatan 2024 + Wang 2025). Para. 2: kontribusi (audit kepala-ke-kepala, E1–E3, temuan utama, pengukuran penuntas). Tulis dalam **Bahasa Indonesia**; tambahkan versi Inggris bila mau portofolio. |
| — | **Kata Kunci** | — | — | 5–7 kata kunci IEEE-style, mis.: *efek Warburg, efisiensi proteom, pemrograman linear, analisis sensitivitas Sobol, metabolisme limpahan, alokasi proteom* |
| **I** | **Pendahuluan** | Pendahuluan (a) | 0,75 | Motivasi Warburg (1 par.); ringkas Shen vs. Kukurugya + Wang (1 par.); **celah yang tersisa**: belum ada audit kepala-ke-kepala + identifiabilitas (1 par.); pernyataan kontribusi + struktur paper (1 par.). Akhiri dengan kalimat kontribusi eksplisit: *"Paper ini memberikan kontribusi: (i)…, (ii)…, (iii)…"* |
| **II** | **Kajian Pustaka** | Pendahuluan (a) — tinjauan pustaka | 0,5 | Tabel atau narasi ringkas: Basan 2015 (hipotesis klasik), Shen 2024, Kukurugya 2024, Wang 2025, sEnz 2024. Untuk tiap sumber: metode, klaim utama, keterbatasan relevan. Akhiri dengan paragraf posisi: "Wang (2025) mendamaikan via model baru; kami mengaudit model asli — pendekatan komplementer." |
| **III** | **Metode** | Metode (b) | 1,25 | **III-A Formulasi Model A** (LP 2-sektor, kendala C1–C5, tujuan — tampilkan persamaan bernomor). **III-B Formulasi Model B** (5-parameter Kukurugya — tampilkan persamaan). **III-C Antarmuka Bersama** (bagaimana keduanya diberi input identik). **III-D Eksperimen E1** (input-bersama). **III-E Eksperimen E2** (batas-akuntansi — definisikan ruang variasi atribusi). **III-F Eksperimen E3** (Sobol/Morris via SALib — definisikan ruang parameter + output). Sertakan **bagan alir** (satu gambar kolom-tunggal). Tambahkan: *"Kode tersedia di [URL repo]."* |
| **IV** | **Hasil** | Hasil & Diskusi (c) — bagian kuantitatif | 1,0 | **IV-A** F1: overlay dua model pada sumbu glukosa bersama — di mana ambang peralihan berbeda, seberapa besar. **IV-B** F2: kurva kapasitas-vs-realisasi — nilai fraksi idle `u_G` tempat verdict membalik (nyatakan sebagai angka eksplisit + pita ketidakpastian). **IV-C** F3: tornado/bar indeks Sobol atas output kontinu — parameter teratas, *"pengukuran penuntas"*. Tiap sub-section: satu gambar + satu kalimat ringkasan kuantitatif. |
| **V** | **Diskusi** | Hasil & Diskusi (c) — bagian kualitatif | 0,6 | Tafsir biologis temuan E1–E3: *mengapa* atribusi enzim penting secara biologi; implikasi "pengukuran penuntas" bagi eksperimentalis. Posisi terhadap Wang: apakah temuan audit mendukung atau mempertanyakan mekanisme heterogenitasnya? Keterbatasan: `ε` dari model mainan, bukan GECKO skala-genom; tiga organisme terbatas. |
| **VI** | **Kesimpulan** | Kesimpulan (d) | 0,25 | 1 paragraf: ringkas temuan (sumber perselisihan = definisi+parameter, bukan struktur; E. coli tidak pernah diperselisihkan; pengukuran X akan menuntaskan). 1 paragraf: kerja masa depan (GECKO/ETFL skala-genom, validasi eksperimental, gandeng heterogenitas Wang). |
| — | **Daftar Pustaka** | Daftar Pustaka (e) | 0,5 | Format **IEEE numbered** `[1] Nama, "Judul," *Jurnal*, vol., hal., tahun.` Minimal 9 referensi sesuai §13. |

---

### Aturan menulis

Tiap persamaan diberi nomor `(1)`, `(2)`, dst. dan dirujuk di teks. Tiap gambar diberi caption dan dirujuk sebagai "Gambar 1", "Gambar 2". Hipotesis dikonfirmasi/dibantah dengan angka eksplisit. **Jangan klaim** "tidak ada yang mendamaikan" — Wang sudah melakukannya; klaim yang benar dan jujur: "belum ada yang mengaudit kedua model asli kepala-ke-kepala." Pertahankan ≤6 halaman ketat — pangkas Kajian Pustaka sebelum memangkas Hasil.

### Pembagian penulisan yang disarankan

| Anggota | Bagian |
|---|---|
| A | §III (Metode) — sudah paling tahu model |
| B | §II (Kajian Pustaka) + §VI (Kesimpulan) + format/referensi |
| C | §IV (Hasil) — sudah pegang data audit |
| D | §I (Pendahuluan) + §V (Diskusi) + gambar final |

Tulis tiap bagian di Markdown dulu, lalu gabungkan. Satu orang bertanggung jawab gabungan + konsistensi notasi.

---

## 11. Luaran 3 — Spesifikasi Dek Slide + Video

Anda ingin ini **terlihat seperti dek** — jadi rancang untuk slide, bukan dinding teks. Target ~12–14 slide, video ≤10 mnt, tiap anggota tampil.

**Sistem desain (konsisten):**
- 16:9. Satu ide per slide. ≤25 kata/slide; gambar yang bicara.
- Palet dua-warna aksen (mis. teal-tua = respirasi, kuning-amber = glikolisis — pakai ulang warna ini di gambar agar dek & laporan selaras visual).
- Satu keluarga sans-serif (Inter / Source Sans). Judul slide besar (32–40 pt), margin lega, tanpa klip-art.
- Tiap slide hasil = satu gambar + satu kalimat-takeaway sebagai judul gaya-pernyataan (mis. *"Verdict membalik hanya dengan mengubah enzim mana yang dihitung"*).

**Slide demi slide:**
1. **Judul** — nama proyek, pemikat ("Dua makalah 2024, jawaban berlawanan — siapa benar?"), tim, mata kuliah.
2. **Teka-teki** — efek Warburg dalam satu diagram (sel memfermentasi walau ada O₂).
3. **Pertarungan 2024 + pendamai 2025** — Shen vs. Kukurugya berdampingan + Wang sebagai "model ketiga"; tabel `ρ`.
4. **Celah yang tersisa** — "Wang menambah model; belum ada yang mengaudit kedua model asli." → RQ + hipotesis.
5. **Model** — gambar dua-sektor: anggaran protein dibagi dua mesin ATP; sebutkan ada dua varian (Shen 2-sektor, Kukurugya 5-parameter).
6. **Tiga eksperimen audit** — E1 input-bersama, E2 batas-akuntansi, E3 identifiabilitas, sebagai ikon.
7. **Bagan alir metode** — diagram arsitektur (§4).
8. **Hasil 1 — overlay dua model (F1)** — di mana mereka benar-benar berbeda.
9. **Hasil 2 — pembalikan batas-akuntansi (F2)** — "sebagian perselisihan itu definisional."
10. **Hasil 3 — pengukuran penuntas (F3, indeks Sobol)** — slide "aha": eksperimen yang menyelesaikan perdebatan.
11. **(Peregangan) cek model-inti (F4)** — peralihan sama di jaringan nyata.
12. **Makna** — sumber perbedaan terlokalisasi; E. coli tak pernah diperselisihkan; pengukuran X menuntaskan; komplementer terhadap Wang.
13. **Keterbatasan & kerja masa depan** — `ε` mainan, skala-genom berikutnya, uji eksperimental.
14. **Penutup / kontribusi** — tautan repo, siapa mengerjakan apa, terima kasih.

**Produksi video:** narasi di atas dek (layar + webcam pojok boleh); sertakan **demo ~60–90 dtk** `python main.py` berjalan dan gambar muncul (spesifikasi mewajibkan demo program singkat); tiap anggota menarasikan satu bagian agar semua tampil; tautkan rekaman di laporan.

---

## 12. Daftar Risiko & Garis-Pangkas

| Risiko | Kemungkinan | Mitigasi / urutan pangkas |
|---|---|---|
| **Premis salah: ternyata kedua makalah "tidak benar-benar bertentangan"** | **Sedang–Tinggi** | Ini justru **temuan, bukan kegagalan.** Bingkai audit sebagai "mendiagnosis apakah & mengapa" sejak awal; bila ternyata semu (kapasitas vs. realisasi), itulah hasil utama yang dapat dipublikasikan. Jangan klaim "kontradiksi" di abstrak — klaim "perselisihan yang dilaporkan". |
| **5 parameter Kukurugya / arah `ρ` salah karena dari parafrase** | Tinggi | **Gerbang Hari-0**: angkat dari tabel/SI asli, verifikasi, tandai E. coli sebagai kasus khusus. Jangan pernah pakai 1,85/0,48/0,32 tanpa verifikasi. |
| **Re-implementasi tidak setia → audit tak sah** | Sedang | **Gerbang validasi (§3.3)**: tiap model wajib reproduksi ≥1 gambar terbit ±20% sebelum dibandingkan. Bila Model B gagal lolos, audit dibatalkan dan paper menjadi laporan replikasi+kegagalan (tetap dapat dinilai). |
| **E2/E3 degenerat (tak ada pembalikan / sensitivitas datar)** | Rendah | Sapu fraksi idle `u_G` cukup lebar; pakai output kontinu (margin/ambang) untuk Sobol, bukan biner; ambil rentang parameter dari SI agar realistis |
| **Sobol salah-pakai pada output biner** | (sudah dikoreksi) | Output E3 **wajib kontinu** (margin `m` atau ambang `g*`); jangan jalankan Sobol atas verdict biner |
| Peregangan M6 (COBRApy) memakan Hari 3–4 | Sedang | **Pangkas M6 lebih dulu.** Inti (Model A+B, E1–E3) sudah proyek lengkap |
| Video lewat 10 mnt | Sedang | Skrip ke ~8 mnt; latih sekali; pangkas slide peregangan |
| Repo tak dapat diakses penilai | Sedang | Uji tautan incognito pada Hari 5 |
| Tinjau pustaka menemukan audit serupa sudah ada | Rendah | Dicek 24 Mei 2026: ada Wang (model baru) & sEnz (metode generik), **tidak ada** re-implementasi kepala-ke-kepala + identifiabilitas atas Shen–Kukurugya. Posisikan jujur terhadap keduanya. |

**Garis-pangkas keras (putuskan di standup Hari-3):** bila F1+F2 belum jadi akhir Hari 3, tinggalkan M6 & semua poles, kirim audit overlay (E1) + batas-akuntansi (E2) saja. Itu sudah memenuhi tiap butir rubrik.

---

## 13. Daftar Pustaka (bawa ke laporan)

1. Shen, Y. et al. *Nat. Chem. Biol.* **20**, 1123–1132 (2024). doi:10.1038/s41589-024-01571-y. (Respirasi lebih efisien-proteom; lindung-nilai proteom.)
2. Kukurugya, M. A., Rosset, S. & Titov, D. V. *PNAS* **121**(46), e2409509121 (2024). doi:10.1073/pnas.2409509121. (Glikolisis lebih cepat per mg protein; model 5-parameter.)
3. Wang, X. *eLife* **13**, RP94586 (2025). doi:10.7554/eLife.94586. (Pendamaian via optimasi pertumbuhan + heterogenitas sel.)
4. *Aerobic glycolysis comes with an enzyme cost but robustness gain.* (News & Views) *Nat. Chem. Biol.* **20**, 1108–1109 (2024). doi:10.1038/s41589-024-01581-w. (Ringkasan independen Shen: per massa enzim, respirasi lebih cepat & lebih efisien-proteom; glikolisis memberi ketahanan hipoksia.)
4. Liu, Y. et al. *sEnz: Sensitivities in protein allocation models.* *Bioinformatics* **40**(12), btae691 (2024). (Kerangka sensitivitas untuk model alokasi-protein via shadow price.)
5. Basan, M. et al. *Nature* **528**, 99–104 (2015). (Metabolisme limpahan / efisiensi proteom.)
6. Pfeiffer, T., Schuster, S. & Bonhoeffer, S. *Science* **292**, 504–507 (2001). (Tukar-tukar hasil-vs-laju ATP.)
7. Herman, J. & Usher, W. *SALib.* *J. Open Source Softw.* **2**(9), 97 (2017). (Indeks Sobol & Morris.)
8. Ebrahim, A. et al. *COBRApy.* *BMC Syst. Biol.* **7**, 74 (2013). (Alat model peregangan.)
9. Salvy, P. & Hatzimanikatis, V. *Nat. Commun.* **11**, 30 (2020). (ETFL — ekstensi skala-genom, kerja masa depan.)

> Cek-silang setiap sitasi terhadap PDF asli sebelum pengumpulan; verifikasi angka `ε` per-organisme di SI tiap makalah, jangan kutip parafrase.

---

## 14. Mengapa Ini Memenangkan Rubrik

| Kriteria penilaian (dari spesifikasi) | Bagaimana proyek ini menjawab |
|---|---|
| **Konsep Biologi** (kedalaman biologi) | Respirasi vs. fermentasi, hasil-vs-laju ATP, Warburg/Crabtree, anggaran proteom, makna lindung-nilai — semua menanggung beban, bukan hiasan. |
| **Analisis Komputasi** (ketelitian) | Dua model bersama eksak (LP + analitik), gerbang validasi terhadap gambar terbit, eksperimen kapasitas-vs-realisasi terkendali, sensitivitas global Sobol+Morris ber-output-kontinu dengan pemeriksaan-silang metode, pipeline reprodusibel + uji. |
| **Inovasi & Problem-Solving** | Bingkai *audit + identifiabilitas* baru (bukan duplikat Wang maupun sEnz), tidak over-engineered — model kecil menjawab pertanyaan meta-ilmiah: *apakah kedua makalah benar-benar bertentangan, dan eksperimen mana yang menuntaskannya.* |
| **Kualitas Laporan** | Struktur dipetakan ke spesifikasi, digerakkan-gambar, kontroversi & posisi terhadap Wang dinyatakan jujur. |
| **Kerja Tim** | Kepemilikan modul jelas + tabel kontribusi. |
| **Presentasi** | Dek judul-pernyataan, gambar-dahulu; tiap anggota tampil; demo langsung. |

---

## 15. Catatan Audit Akademik (changelog kritik — disimpan demi transparansi)

Bagian ini didokumentasikan agar keputusan desain dapat dipertanggungjawabkan. Spec ini telah dikritik-sendiri dan direvisi pada 24 Mei 2026 setelah verifikasi sumber primer (Shen 2024; Kukurugya 2024; companion *Nat. Chem. Biol.* 1108; Wang 2025; sEnz 2024). Tujuh celah ditemukan dan diperbaiki:

1. **Risiko premis (kritis).** Versi awal membingkai Shen↔Kukurugya sebagai kontradiksi langsung. Verifikasi menunjukkan keduanya kemungkinan **mengukur besaran berbeda** — kapasitas maksimal `V·γ` (Kukurugya) vs. fluks ter-realisasi `J/m_enzim` (Shen). *Perbaikan:* RQ1 & Hipotesis ditulis ulang; "kontradiksi" diganti "perselisihan yang dilaporkan"; ini dijadikan hipotesis utama, bukan asumsi.
2. **Model A tidak setia (kritis).** Versi awal memakai kapasitas **O₂** sebagai kendala pengikat — keliru, karena Warburg terjadi saat O₂ ada. *Perbaikan:* Model A dibangun ulang dengan **ketersediaan glukosa + anggaran proteom** sebagai knob, selaras dengan Kukurugya, sehingga E1 punya sumbu bersama yang sah.
3. **Presisi palsu (integritas).** Tabel `ρ` (1,85/0,48/0,32) diturunkan dari parafrase, berpotensi salah tanda untuk E. coli (kasus Pta-AckA). *Perbaikan:* ditandai provisional + **gerbang Hari-0** mewajibkan pengambilan dari tabel/SI asli.
4. **E2 salah sasaran.** "Batas-akuntansi enzim" diturunkan dari status utama. *Perbaikan:* E2 utama kini **kapasitas vs. realisasi** (sapu fraksi idle `u_G`), sesuai temuan ekspresi-berlebih glikolitik Shen; akuntansi-enzim jadi sub-eksperimen sekunder.
5. **E3 tidak well-posed.** Sobol dijalankan atas output biner ("verdict"). *Perbaikan:* output diubah jadi **kontinu** (margin efisiensi `m` atau ambang glukosa `g*`); distribusi input diambil dari rentang SI.
6. **Tak ada gerbang validasi.** *Perbaikan:* ditambahkan **gerbang validasi wajib** — tiap model harus reproduksi ≥1 gambar terbit ±20% sebelum masuk audit.
7. **Kabar baik terlewat.** Model Kukurugya ternyata **analitik, solusi unik, tanpa parameter bebas** (5 parameter). *Perbaikan:* Model B didefinisikan presisi; beban implementasi turun, tapi gerbang konsep Hari-0 ditambah.

**Sikap akademik yang wajib dijaga.** (a) Wakili kedua model **secara murah hati** (charitable) — pakai definisi masing-masing pihak; audit **mendiagnosis sumber perbedaan**, bukan menyatakan siapa "salah". (b) Setiap klaim kuantitatif (mis. titik-balik `u_G`, ambang `g*`) **wajib disertai pita ketidakpastian** dari rentang parameter, bukan angka tunggal. (c) Nyatakan keterbatasan jujur: parameter dari model coarse-grained/SI, bukan pengukuran sendiri; tiga organisme; in-silico, bukan validasi basah.

**Risiko tersisa yang tak bisa dihilangkan (akui di paper §V):** bila parameter SI tak lengkap untuk salah satu organisme, batasi audit ke organisme yang datanya lengkap dan nyatakan demikian — lebih baik sempit-tapi-sah daripada luas-tapi-menebak.