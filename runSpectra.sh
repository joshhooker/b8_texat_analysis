root -l $1 << EOF
.X MakeChain.C
.L Spectra.C+
Spectra t(mfmData)
t.Loop()
EOF