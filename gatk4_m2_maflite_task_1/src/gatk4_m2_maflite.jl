#!/usr/local/bin/julia
# ARGS=["REBC-AC8L-TP","REBC-AC8L-NT","tmp1.tsv","REBC-AC8L-TP-NT.M2_maflite.tsv"]
# ARGS=["SU2CLC-MGH-1047-TM-01","SU2CLC-MGH-1047-BL-01","tmp3.tsv","RP-1066_SU2CLC-MGH-1047"]
# ARGS=["753TD","753ND","tmp3.tsv","753TN"]

using DataFrames
using DelimitedFiles
using Statistics

tumor_id=ARGS[1]
normal_id=ARGS[2]
file1=ARGS[3]
pair_id=ARGS[4]
file1a=pair_id*".raw.tsv"

isString(x::Number)=false
#isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

df = readtable(file1)




# size(df)
# describe(df)
# print(df)
FIELDS=names(df)
FIELDX=[:x,:QUAL,:ID,:NORMAL_FT,:NORMAL_MMQ,:NORMAL_PID, :NORMAL_MBQ, :NORMAL_GQ, :NORMAL_OBF, :NORMAL_OBQRC,:NORMAL_AF,:NORMAL_FOXOG,
 :NORMAL_SA_MAP_AF, :NORMAL_MPOS, :NORMAL_PL, :NORMAL_OBAM, :NORMAL_OBQ, :NORMAL_ALT_F2R1, :NORMAL_SA_POST_PROB, :NORMAL_REF_F2R1,
 :NORMAL_GT, :NORMAL_MCL, :NORMAL_MFRL, :NORMAL_OBP, :NORMAL_OBAMRC, :NORMAL_PGT, :NORMAL_REF_F1R2,  :NORMAL_ALT_F1R2,
 :TUMOR_FT, :TUMOR_PID,  :TUMOR_GQ, :TUMOR_OBF, :TUMOR_OBQRC, :TUMOR_SA_MAP_AF, :TUMOR_PL, :TUMOR_OBAM,
 :TUMOR_OBQ, :TUMOR_SA_POST_PROB,:TUMOR_GT,:TUMOR_OBP, :TUMOR_OBAMRCM, :TUMOR_PGT]

for i=1:length(FIELDX)
    if (FIELDX[i] in FIELDS)
        delete!(df, FIELDX[i])
    end
end
#delete!(df, [:FILTER,:QUAL, :NORMAL_PID,:x,:TUMOR_GT,:TUMOR_PL,:TUMOR_GQ,:TUMOR_DP,:TUMOR_PGT,:NORMAL_FOXOG,:NORMAL_QSS,:NORMAL_ALT_F2R1,:NORMAL_ALT_F1R2,:NORMAL_REF_F2R1,:NORMAL_REF_F1R2,:NORMAL_GT,:NORMAL_PGT,:NORMAL_GQ,:NORMAL_PL,:NORMAL_DP])
#if ( :TUMOR_FOXOG in names(df) )  # FOXOG messed up by VariantAnnotator
#    delete!(df, [:TUMOR_FOXOG])
#end
# print(df)


for c in names(df)
    println(c)
    k=ismissing(df[c])
    if Statistics.mean(k)==1.0
        delete!(df, c)
        continue
    end
    # if isa(df[c],DataArray{Float64,1})  & (mean(k)>0)
    #     v=df[c]
    #     v=map(x -> replace(x,NA,NaN), v)
    #     df[c]=v
    # end
    if isa(df[c],Array{String,1}) 
        df[k,c] = ""
    end
end
#rename!(df, [:normal_tumor_alt_count, :normal_tumor_ref_count], [:t_alt_count, :t_ref_count])
head(df)
a= df[:CHRO]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]"=> "23"), a)
    a=map(x -> replace(x,r"[Y]"=> "24"), a)
    a=map(x -> replace(x,r"[MT]"=> "25"), a)
    a=map(x -> replace(x,r"[M]"=> "25"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
# print(df)
sort!(df, cols = [:a, :POS])
delete!(df, [:a])

println("")
println(file1a)
println("")

for c in [:TUMOR_F1R2_ALT,:TUMOR_F2R1_ALT,:TUMOR_F1R2_REF,:TUMOR_F2R1_REF]
    if ~isString(df[1,c])
        df[c] = map(x -> string(x),df[c])
    end
    k=findall(map(x->x=="NA",df[c]))
    df[k,c] = ""
end


open(file1a, "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Array,df), '\t')
end

a=df[:ALT]
r=df[:REF]
n=length(a)
ar=fill("",n)
for i in 1:n
    #println(i)
    ar[i]=string(a[i],":",r[i])
end

l1=map(x -> length(x), ar)
df[:INDEL]=map(x -> x>3,l1)
df[:Variant_Type]=fill("SNV",n)
df[df[:INDEL],:Variant_Type]="INDEL"
dfindel=df[df[:INDEL],:]
dfsnv=df[.~df[:INDEL],:]
#delete!(dfsnv, [:INDEL])
#delete!(dfindel, [:INDEL])
#delete!(df, [:INDEL])
df2=df

# set ALT0 to original ALT,shorten ALT to at most 5 bases, then copy back after merge

df2[:ALT0]=df2[:ALT]
alt0=df2[:ALT0]
alt=fill("",n)
n=length(alt0)
l0=map(x -> length(x), alt0)
for i in 1:n
    l1=min(l0[i],5)
    alt1=string(alt0[i])
    alt[i]=alt1[1:l1]
end
df2[:ALT]=alt



delete!(df, [:INDEL, :Variant_Type,:ALT0])



# label M1 with M1 flag
a= df[:CHRO]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]"=> "23"), a)
    a=map(x -> replace(x,r"[Y]"=> "24"), a)
    a=map(x -> replace(x,r"[MT]"=> "25"), a)
    a=map(x -> replace(x,r"[M]"=> "25"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
#print(df)
sort!(df, cols = [:a, :POS])
delete!(df, [:a])



for c in names(df)
    println(c)
    #k=findall(map(x->ismissing(x),df[c]))
    k1=map(x->ismissing(x)*1,df[c])
    k=findall(k1.>0)
    k2=findall(k1.==0)
        #k=fomismissing(df[c])
    if Statistics.mean(k1)==1.0
        delete!(df, c)
        continue
    end
    println(typeof(df[k2[1],c]))
    if isString(df[k2[1],c])
        df[k,c] = ""
    end
end

# maflite fields 
# SNP: build    chr start   end ref_allele  tum_allele1 tum_allele2 tumor_barcode   normal_barcode  tumor_f init_t_lod  t_lod_fstar t_alt_count t_ref_count judgement
# INDEL: build  contig  start_position  end_position    ref_allele  alt_allele1 alt_allele2 tumor_name  normal_name tumor_f n_ref_count n_alt_count t_ref_count t_alt_count init_n_lod  init_t_lod  judgement

df[:end]=df[:POS]

# start for indels is POS+1
kdel = map(x-> length(x)>1, df[:REF])
p=df[:POS]
p=p+kdel
df[:POS]=p
k = findall(kdel)
ref=df[k,:REF]
ref=map(x->x[2:end],ref)
alt=df[k,:ALT]
alt=map(x->x[2:end],alt)
k1 = findall(map(x-> length(x)<1, alt))
alt[k1]=map(x-> "-", alt[k1])
df[k,:REF]=ref
df[k,:ALT]=alt
dref = map(x-> length(x), ref)
dalt = map(x-> length(x), alt)
ddel = dref
df[k,:end]=df[k,:POS]+ddel




kins = map(x-> length(x)>1, df[:ALT])
# p=df[:POS]
# p=p+kins
# df[:POS]=p
# df[find(kins),:POS]
k = findall(kins)
ref=df[k,:REF]
ref=map(x->x[2:end],ref)
k1 = findall(map(x-> length(x)<1, ref))
ref[k1]=map(x-> "-", ref[k1])
alt=df[k,:ALT]
alt=map(x->x[2:end],alt)
df[k,:REF]=ref
df[k,:ALT]=alt
df[k,:end]=df[k,:POS].+1



# add build column
df[:build]=fill("37",size(df,1))
df[:tumor_barcode]=fill(tumor_id,size(df,1))
df[:normal_barcode]=fill(normal_id,size(df,1))
df[:judgement]=fill("KEEP",size(df,1))


# M1 counts
df[:n_alt_count]=df[:NORMAL_AD_ALT]
df[:n_ref_count]=df[:NORMAL_AD_REF]
df[:t_lod_fstar]=df[:TLOD]
df[:tumor_f]=df[:TUMOR_AF]
df[:t_alt_count]=df[:TUMOR_AD_ALT]
df[:t_ref_count]=df[:TUMOR_AD_REF]
rename!(df, [:CHRO,:POS, :REF, :ALT], [:chr, :start, :ref_allele,:alt_allele])#from = [:DB, :HCNT, :MAX_ED,:RPA,:STR, :RU, :NLOD,:TUMOR_ALT_F1R2,:MIN_ED,:ECNT, :TLOD, :TUMOR_PID,:TUMOR_AF,:TUMOR_FOXOG, :TUMOR_QSS, :TUMOR_ALT_F2R1,:TUMOR_AD_REF,:TUMOR_AD_ALT, :TUMOR_REF_F2R1, :TUMOR_REF_F1R2,:NORMAL_AF,:NORMAL_AD_REF,:NORMAL_AD_ALT]
#from = ["DB","HCNT","MAX_ED","RPA","STR","RU","NLOD","TUMOR_ALT_F1R2","MIN_ED","ECNT","TLOD","TUMOR_PID","TUMOR_AF","TUMOR_FOXOG","TUMOR_QSS","TUMOR_ALT_F2R1","TUMOR_AD_REF","TUMOR_AD_ALT","TUMOR_REF_F2R1","TUMOR_REF_F1R2","NORMAL_AF","NORMAL_AD_REF","NORMAL_AD_ALT"]
#println(from)

#for c in names(df)
#    c1 = string(c)
#    println(c1)
#    if length(find(Bool[contains(c1,i) for i in from]))>0
#        m2c = Symbol(string("M2_",c1))
#        rename!(df,c,m2c)
#    end
#end

maf=df
for c in names(maf)
    #println(c)
    if ~isString(maf[1,c])
        maf[c] = map(x -> string(x),maf[c])
    end
end


for c in names(maf)
    println(c)
    k1=map(x->ismissing(x)*1,maf[c])
    k=findall(k1.>0)
    k2=findall(k1.==0)
    if mean(k)==1.0
        delete!(maf, c)
        continue
    end
    k=findall(map(x->x=="NA",maf[c]))
    maf[k,c] = ""
    k=findall(map(x->uppercase(x)=="TRUE",maf[c]))
    maf[k,c] = "1"
    k=findall(map(x->uppercase(x)=="FALSE",maf[c]))
    maf[k,c] = "0"
end

maflite_all=pair_id*".m2.all.maflite.tsv"
open(maflite_all, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end

kpass = findall(map(x-> uppercase(x)=="PASS", maf[:FILTER]))
maf=maf[kpass,:]
delete!(maf,[:FILTER])


maflite_pass=pair_id*".m2.pass.maflite.tsv"
open(maflite_pass, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end

