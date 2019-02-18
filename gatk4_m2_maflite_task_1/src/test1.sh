cd /opt/gatk4_m2_maflite/gatk4_m2_maflite_task_1/src
TID="753TD"
NID="753ND"
PID="753TN"
VCF1="/opt/test/753_pair.MuTect2.call_stats.vcf"
#export JULIA_LOAD_PATH="."
julia vcf2txt.jl $VCF1  tmp1.tsv

sed "s/${TID}/TUMOR/g" tmp1.tsv > tmp2.tsv
sed "s/${NID}/NORMAL/g" tmp2.tsv > tmp3.tsv

julia gatk4_m2_maflite.jl $TID $NID tmp3.tsv $PID


