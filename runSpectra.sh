root -l $1 << EOF
.L Spectra.C+
Spectra t(mfmData)
t.Loop()
EOF