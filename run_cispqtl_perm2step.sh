#!/bin/bash
#---------------------------------------#---------------------------------------
dir_wk=.;
[ -d ${dir_wk}/output ] || mkdir -p ${dir_wk}/output;
[ -d ${dir_wk}/output/out_perm1k ] || mkdir -p ${dir_wk}/output/out_perm1k;
[ -d ${dir_wk}/output/out_perm10k ] || mkdir -p ${dir_wk}/output/out_perm10k;
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     plink_=${dir_wk}/plink.v1.90b.linux;;
    Darwin*)    plink_=${dir_wk}/plink.v1.90b.macos;;
    *)          exit
esac;
threads_=8; # For hyperthreading in plink
bpdist_=1000000;
#-Genotype, protein abundance, covariate data and protein position (from gencode v26)
file_geno=${dir_wk}/data/example_data_geno;
file_prt=${dir_wk}/data/example_data_prt.txt;
file_covar=${dir_wk}/data/example_data_covar.txt;
file_gencode=${dir_wk}/data/example_data_gencode.txt;
ncovar2_=`awk -F" " 'NR == 1{print NF; exit;}' ${file_covar}`;
#-Parameters of permutations for first round (1k)
mperm1_=1000;
seed1_=2024120914;
file1k=${dir_wk}/output/out_perm1k_summ.txt;
#-Parameters of permutations for second round (10k)
thr_cnt1k=25;
file10k_input=${dir_wk}/output/input_perm10k.txt;
file10k=${dir_wk}/output/out_perm10k_summ.txt;
mperm2_=10000;
seed2_=2412091812;
#-Output cis-pQTL results
fileout=${dir_wk}/output/out_cispqtl.csv;
#---------------------------------------#---------------------------------------
#-Step 1. 1000 permutations for cis-pQTL analysis
start0=1;
end0=$(wc -l < ${file_gencode});
for ((kk = ${start0}; kk <= ${end0}; kk += 1)); do
	genepos0=($(sed -n "${kk}p" ${file_gencode}));
	#-Autosome and sample size no less than number of covariates
	if [[ ${genepos0[1]} =~ ^[0-9]+$ ]] && [ ${genepos0[5]} -gt ${ncovar2_} ]; then
		s0a=$((genepos0[2] - bpdist_));
		s0a=$((s0a > 0 ? s0a : 0));
		fileout0_=${dir_wk}/output/out_perm1k/${genepos0[0]};
		if [ ! -f "${fileout0_}.minpval.txt" ] || [ $(wc -l < "${fileout0_}.minpval.txt") -lt $(( mperm1_ + 1 )) ]; then
			${plink_} --threads ${threads_} --seed ${seed1_} --bfile ${file_geno} --chr ${genepos0[1]} --from-bp ${s0a} --to-bp $((genepos0[3] + bpdist_)) --allow-no-sex --linear hide-covar mperm=${mperm1_} --mperm-save-all --covar ${file_covar} --pheno ${file_prt} --pheno-name ${genepos0[0]} --out ${fileout0_} > /dev/null;
			sed -i 's/ \+/ /g; s/^ //g; s/ $//g;' ${fileout0_}.assoc.linear;
			sed 's/ \+/ /g; s/^ //g; s/ $//g;' ${fileout0_}.assoc.linear.mperm | cut -d" " -f 3 | paste -d " " ${fileout0_}.assoc.linear - > ${fileout0_}.assoc.linear.emp1;
#----
R --quiet --vanilla << Rcmd
	require(data.table, quiet = T);
	ncnt_R <- as.integer("${ncovar2_}");
	out0_R <- "${fileout0_}";
	nmiss0 <- fread(paste0(out0_R, ".assoc.linear.emp1"), head = T, sep = " ", select = 6)[[1]] - ncnt_R;
	id0 <- nmiss0 > 0;
	nmiss0 <- nmiss0[id0];	
	dump0 <- fread(paste0(out0_R, ".mperm.dump.all"), head = F, drop = 1)[, id0, with = F];
	grp0 <- unique(nmiss0);
	pval0 <- sapply(grp0, function(x0) {
		dump2 <- dump0[, nmiss0 %in% x0, with = F];
		stat0 <- apply(dump2, 1, max, na.rm = T);
		return(2 * pt(-stat0, df = x0));
	});
	pval2 <- signif(apply(pval0, 1, min, na.rm = T), 6);
	write(pval2, file = paste0(out0_R, ".minpval.txt"), ncolumns = 1);
Rcmd
#----
			rm ${fileout0_}.{log,nosex,mperm.dump.all,assoc.linear,assoc.linear.mperm};
		fi;
		cnt1k=`awk -v nr0=$((mperm1_ + 1)) 'BEGIN{cnt = 0;} NR == 1{tt = $0; next;} NR >= 2 && $0 <= tt {cnt++;} END{if(NR != nr0 || tt == 2) cnt = "NA"; print cnt;}' ${fileout0_}.minpval.txt`;
		snp1k=`awk -v start=${genepos0[2]} -v end=${genepos0[3]} -F" " 'function dist0(x0f, xminf, xmaxf) {d0f=0; if(x0f<xminf) d0f=xminf-x0f; if(x0f>xmaxf) d0f=x0f-xmaxf; return d0f;} BEGIN{pmin = 1;} NR >= 2 && $9 != "NA"{if($9 < pmin) {out = $0; pmin = $9; cnt = 1; d0 = dist0($3, start, end); next;} if($9 == pmin) {cnt++; d1 = dist0($3, start, end); if(d1 < d0) {out = $0; d0 = d1;}}} END{print out" "cnt" "d0;}' ${fileout0_}.assoc.linear.emp1`;
		echo "OUTPUT_NEEDED: "${genepos0[@]}" "${cnt1k}" "${snp1k};
	fi;
done > ${dir_wk}/output/__tmp1k__.txt;
[ -f ${file1k} ] | grep "^OUTPUT_NEEDED: " ${dir_wk}/output/__tmp1k__.txt | sed 's/.*OUTPUT_NEEDED: //g' > ${file1k};
#---------------------------------------#---------------------------------------
#-Step 2. 10000 permutations for cis-pQTL analysis
# Keep proteins with p-value <= (25+1)/(1000+1) in 1k permutations
cut -d" " -f1-8 ${file1k} | awk -F" " -v thr=${thr_cnt1k} '$7 <= thr {print $0;}' > ${file10k_input};
file_gencode=${file10k_input};
start0=1;
end0=$(wc -l < ${file_gencode});
for ((kk = ${start0}; kk <= ${end0}; kk += 1)); do
	genepos0=($(sed -n "${kk}p" ${file_gencode}));
	if [[ ${genepos0[1]} =~ ^[0-9]+$ ]] && [ ${genepos0[5]} -gt ${ncovar2_} ]; then
		s0a=$((genepos0[2] - bpdist_));
		s0a=$((s0a > 0 ? s0a : 0));
		fileout0_=${dir_wk}/output/out_perm10k/${genepos0[0]};
		if [ ! -f "${fileout0_}.minpval.txt" ] || [ $(wc -l < "${fileout0_}.minpval.txt") -lt $(( mperm2_ + 1 )) ]; then
			${plink_} --threads ${threads_} --seed ${seed2_} --bfile ${file_geno} --chr ${genepos0[1]} --from-bp ${s0a} --to-bp $((genepos0[3] + bpdist_)) --allow-no-sex --linear hide-covar mperm=${mperm2_} --mperm-save-all --covar ${file_covar} --pheno ${file_prt} --pheno-name ${genepos0[0]} --out ${fileout0_} > /dev/null;
			sed -i 's/ \+/ /g; s/^ //g; s/ $//g;' ${fileout0_}.assoc.linear;
			sed 's/ \+/ /g; s/^ //g; s/ $//g;' ${fileout0_}.assoc.linear.mperm | cut -d" " -f 3 | paste -d " " ${fileout0_}.assoc.linear - > ${fileout0_}.assoc.linear.emp1;
#----
R --quiet --vanilla << Rcmd
	require(data.table, quiet = T);
	ncnt_R <- as.integer("${ncovar2_}");
	out0_R <- "${fileout0_}";
	nmiss0 <- fread(paste0(out0_R, ".assoc.linear.emp1"), head = T, sep = " ", select = 6)[[1]] - ncnt_R;
	id0 <- nmiss0 > 0;
	nmiss0 <- nmiss0[id0];	
	dump0 <- fread(paste0(out0_R, ".mperm.dump.all"), head = F, drop = 1)[, id0, with = F];
	grp0 <- unique(nmiss0);
	pval0 <- sapply(grp0, function(x0) {
		dump2 <- dump0[, nmiss0 %in% x0, with = F];
		stat0 <- apply(dump2, 1, max, na.rm = T);
		return(2 * pt(-stat0, df = x0));
	});
	pval2 <- signif(apply(pval0, 1, min, na.rm = T), 6);
	write(pval2, file = paste0(out0_R, ".minpval.txt"), ncolumns = 1);
Rcmd
#----
			rm ${fileout0_}.{log,nosex,mperm.dump.all,assoc.linear,assoc.linear.mperm};
		fi;
		cnt10k=`awk -v nr0=$((mperm2_ + 1)) 'BEGIN{cnt = 0;} NR == 1{tt = $0; next;} NR >= 2 && $0 <= tt {cnt++;} END{if(NR != nr0 || tt == 2) cnt = "NA"; print cnt;}' ${fileout0_}.minpval.txt`;
		snp10k=`awk -v start=${genepos0[2]} -v end=${genepos0[3]} -F" " 'function dist0(x0f, xminf, xmaxf) {d0f=0; if(x0f<xminf) d0f=xminf-x0f; if(x0f>xmaxf) d0f=x0f-xmaxf; return d0f;} BEGIN{pmin = 1;} NR >= 2 && $9 != "NA"{if($9 < pmin) {out = $0; pmin = $9; cnt = 1; d0 = dist0($3, start, end); next;} if($9 == pmin) {cnt++; d1 = dist0($3, start, end); if(d1 < d0) {out = $0; d0 = d1;}}} END{print out" "cnt" "d0;}' ${fileout0_}.assoc.linear.emp1`;
		echo "OUTPUT_NEEDED: "${genepos0[@]}" "${cnt10k}" "${snp10k};
	fi;
done > ${dir_wk}/output/__tmp10k__.txt;
[ -f ${file10k} ] | grep "^OUTPUT_NEEDED: " ${dir_wk}/output/__tmp10k__.txt | sed 's/.*OUTPUT_NEEDED: //g' > ${file10k};
#---------------------------------------#---------------------------------------
#-Step 3. Merge 1k permutation and 10k permutation and calculate q-values
cut -d" " -f 1-7,10-12,14-17 ${file1k} | awk -F" " '$7 != "NA" {print $0;}' | awk -v nperm1=${mperm1_} -v nperm2=${mperm2_} -F" " 'BEGIN{print "gene_id chr start end strand nobs in11k SNP_ID SNP_BP A1 NMISS BETA STAT P pvalue"} NR==FNR{arr[$1]=(1 + $7 + $8)/(1 + nperm1 + nperm2); next;} {pval=(1 + $7)/(1 + nperm1); $7 = 0; if($1 in arr) {pval = arr[$1]; $7=1;} printf $0; printf " %0.8g\n", pval;}' ${file10k} - > ${dir_wk}/output/__tmp11k__.txt;
R --quiet --vanilla << Rcmd
	x0a <- sapply(c("data.table", "qvalue"), require, char = T, quiet = T);
	x0b <- fread("${dir_wk}/output/__tmp11k__.txt", head = T);
	x0c <- cbind(x0b, qvalue = qvalue(x0b[["pvalue"]])[["qvalues"]]);
	write.table(x0c, file = "${fileout}", sep = ",", quote = F, row.names = F, col.names = T);
Rcmd
#---------------------------------------#---------------------------------------
#-
rm ${dir_wk}/output/__tmp{1,10,11}k__.txt;
#-