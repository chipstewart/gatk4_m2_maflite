#!/usr/local/bin/julia
# ARGS=["753TD","753ND","tmp3.tsv","753_pair"]

using DataFrames
using CSV
using DelimitedFiles
using Statistics

tumor_id=ARGS[1]
normal_id=ARGS[2]
file1=ARGS[3]
pair_id=ARGS[4]
file1a=pair_id*".raw.tsv"

isString(x::Number)=false
isString(x::Missing)=false
isString(x::AbstractString)=true

df=CSV.File(file1,delim='\t') |> DataFrame


# size(df)
# describe(df)
# print(df)
FIELDS=names(df)
FIELDSK=[:CHRO, :POS, :REF, :ALT, :FILTER, :POPAF, :NLOD, :PON, :TLOD, :NORMAL_AD_REF, :NORMAL_AD_ALT, :TUMOR_AD_REF, :TUMOR_AD_ALT]

for i=1:length(FIELDS)
    if ~(FIELDS[i] in FIELDSK)
        deletecols!(df, FIELDS[i])
    end
end

df[:VCF_REF]=map(x->x,df[:REF])
df[:VCF_ALT]=map(x->x,df[:ALT])
df[:VCF_POS]=map(x->x,df[:POS])

for c in names(df)
    println(c)
    k=ismissing(df[c])
    if k == false
        continue
    end
    if Statistics.mean(k)==1.0
        deletecols!(df, c)
        continue
    end
    if isa(df[c],Array{String,1}) 
        df[k,c] = ""
    end
end
a= df[:CHRO]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]"=> "23"), a)
    a=map(x -> replace(x,r"[Y]"=> "24"), a)
    a=map(x -> replace(x,r"[MT]"=> "25"), a)
    a=map(x -> replace(x,r"[M]"=> "25"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
sort!(df, [:a, :POS])
deletecols!(df, [:a])

println("")
println(file1a)
println("")

open(file1a, "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Matrix,df), '\t')
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



deletecols!(df, [:INDEL, :Variant_Type,:ALT0])



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
sort!(df, [:a, :POS])
deletecols!(df, [:a])



for c in names(df)
    println(c)
    #k=findall(map(x->ismissing(x),df[c]))
    k1=map(x->ismissing(x)*1,df[c])
    k=findall(k1.>0)
    k2=findall(k1.==0)
        #k=fomismissing(df[c])
    if Statistics.mean(k1)==1.0
        deletecols!(df, c)
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
# multiple alt alleles
df[:multiALT]=map(x->1*(length(split(x,":"))>1),df[:ALT])
k=findall(df[:multiALT].>0)
df[:ALT]=map(x->split(x,":")[1],df[:ALT])
print(join(names(df),"\n"))
df[:POPAF]=map(x->split(x,":")[1],df[:POPAF])
df[:NLOD]=map(x->split(x,":")[1],df[:NLOD])
df[:TLOD]=map(x->split(x,":")[1],df[:TLOD])

# MNPs
df[:MNP]=1 .*[(length(r)>1)&(length(r)==length(a))  for (r,a) = zip(df[:REF], df[:ALT])]
LA=map(x->length(x),df[:ALT])
LR=map(x->length(x),df[:REF])
k=findall(df[:MNP].>0)
df[k,:end]=df[k,:POS]+LA[k].-1

# start for indels is POS+1
df[:DEL]=1 .*[(length(r)>1)&(length(r)>length(a))  for (r,a) = zip(df[:REF], df[:ALT])]
p=df[:POS]
p=p+df[:DEL]
df[:POS]=p
k = findall(df[:DEL].>0)
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



# calc end for INS
df[:INS]=1 .*[(length(a)>1)&(length(a)>length(r))  for (r,a) = zip(df[:REF], df[:ALT])]
k = findall(df[:INS].>0)
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
df[:t_alt_count]=df[:TUMOR_AD_ALT]
df[:t_ref_count]=df[:TUMOR_AD_REF]
#rename!(df, [:CHRO,:POS, :REF, :ALT], [:chr, :start, :ref_allele,:alt_allele])#from = [:DB, :HCNT, :MAX_ED,:RPA,:STR, :RU, :NLOD,:TUMOR_ALT_F1R2,:MIN_ED,:ECNT, :TLOD, :TUMOR_PID,:TUMOR_AF,:TUMOR_FOXOG, :TUMOR_QSS, :TUMOR_ALT_F2R1,:TUMOR_AD_REF,:TUMOR_AD_ALT, :TUMOR_REF_F2R1, :TUMOR_REF_F1R2,:NORMAL_AF,:NORMAL_AD_REF,:NORMAL_AD_ALT]
rename!(df, [f => t  for (f,t) = zip([:CHRO,:POS, :REF, :ALT], [:chr, :start, :ref_allele,:alt_allele])])

maf=df

permutecols!(maf, [:chr,:start,:end,:ref_allele,:alt_allele,:tumor_barcode,:normal_barcode,:POPAF,:NLOD,:PON,:TLOD,:build,:judgement,:n_alt_count,:n_ref_count,:t_lod_fstar,:t_alt_count,:t_ref_count,:NORMAL_AD_REF,:NORMAL_AD_ALT,:TUMOR_AD_REF,:TUMOR_AD_ALT,:FILTER,:multiALT, :MNP, :DEL, :INS,:VCF_REF,:VCF_ALT,:VCF_POS])

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
    if mean(k1)==1.0
        deletecols!(maf, c)
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
    #writedlm(f, convert(Array,maf), '\t')
    writedlm(f, convert(Matrix,maf), '\t')
end

kpass = findall(map(x-> uppercase(x)=="PASS", maf[:FILTER]))
maf=maf[kpass,:]
deletecols!(maf,[:FILTER])


maflite_pass=pair_id*".m2.pass.maflite.tsv"
open(maflite_pass, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Matrix,maf), '\t')
end

