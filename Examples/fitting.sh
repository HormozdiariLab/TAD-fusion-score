### fit 40kb for Schmitt
for i in `seq 22 22`
do
  chr="$i"
  if [ "$i" = "23" ]
  then
    chr="X"
  fi
  resolution=40000
  min_length=1
  mid_length=-1
  max_length=10
  seg_length=600
  min_seg_overlap=50
  zero_val=0.1
  debug=0
  ./../src/fit_hic_model -fn ../Data/GM12878.40Kb.raw.chr22.mat -ff full_matrix_format -res ${resolution} -mn ${min_length} -md ${mid_length} -mx ${max_length} -method segmentation -sg ${seg_length} -mso ${min_seg_overlap} -zero ${zero_val} -debug ${debug} -cm "chr${chr}" -of Output/GM12878.40kb.chr${i}.model >> Output/fit_40kb_history
done


