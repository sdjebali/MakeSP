(* bio_objects.ml *)



type strand = | Forward  
              | Reverse



module Exon =
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	score: float;
	cat: string;
	trid: string;
	gnid: string;
  }

  let create chr gb ge st sc ct tr gn = { chrom=chr; gbeg=gb; gend=ge; str=st; score= sc; cat=ct; trid=tr; gnid=gn}
  let null = create "" (-1) (-1) Forward  0.0 "" "" ""

  let compare e1 e2 = 
    if ((Pervasives.compare e1.gnid e2.gnid) !=0) then
      Pervasives.compare e1.gnid e2.gnid
    else
      begin
	if ((Pervasives.compare e1.gbeg e2.gbeg) !=0) then
	  Pervasives.compare e1.gbeg e2.gbeg
	else
	  Pervasives.compare e1.gend e2.gend
      end
    
  let chrom e = e.chrom
  let gbeg e = e.gbeg
  let gend e = e.gend
  let str e = e.str
  let score e = e.score
  let cat e = e.cat
  let trid e = e.trid
  let gnid e = e.gnid
end 


module Transcript =
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	cat: string;
	exlist: Exon.t list;
	trid: string;
	gnid: string;
  }
	  
  let create chr gb ge st ct el tr gn = { chrom=chr; gbeg=gb; gend=ge; str=st; cat=ct;  exlist=el; trid=tr; gnid=gn}

  let chrom t = t.chrom
  let gbeg t = t.gbeg
  let gend t = t.gend
  let str t = t.str	
  let cat t = t.cat
  let exlist t = t.exlist
  let trid t = t.trid
  let gnid t = t.gnid
end



module Gene = 
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	trarr: Transcript.t array;
	exarr: Exon.t array;
	gnid: string;
  }

  let create chr gb ge st ta ea gn = { chrom=chr; gbeg=gb; gend=ge; str=st; trarr=ta; exarr= ea; gnid=gn}
	  
  let chrom g = g.chrom
  let gbeg g = g.gbeg
  let gend g = g.gend
  let str g = g.str	
  let trarr g = g.trarr
  let exarr g = g.exarr
  let gnid g = g.gnid
end



(* Exon projection. In a gene it is a maximal set of overlapping exons *)
module ExonProj =
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	nbex: int;  (* number of exons if is made of *)
	exarr: Exon.t array;  (* The exons it is made of *)
	noingn: int;    (* number of the exon projection among all the exon projections of the gene (from 5') *)
	begingn: int;   (* begining of the exon projection in the virtual cDNA made by joining all the exon projections of the gene *)
	endingn: int;   (* end of the exon projection in the virtual cDNA made by joining all the exon projections of the gene *)
 }
	  
  let create chr gb ge st ne ea no bgn egn = { chrom=chr; gbeg=gb; gend=ge; str=st; nbex=ne; exarr=ea; noingn=no; begingn=bgn; endingn=egn}
  
  let chrom e = e.chrom
  let gbeg e = e.gbeg
  let gend e = e.gend
  let str e = e.str
  let nbex e = e.nbex
  let exarr e = e.exarr
  let noingn e = e.noingn
  let begingn e = e.begingn 
  let endingn e = e.endingn
end 




(* Segmented (Exon) projection. The different exons of an exon projection define elementary segments called Segmented (Exon) projections *)
module SegProj = 
struct
  type t = {
	chrom: string;
	gbeg: int;
	gend: int;
	str: strand;
	score : float;
	exproj: ExonProj.t;
	coverage: int;   (* number of transcripts it belongs to *)
	noinexproj: int; (* number of the segmented projection in the exon projection from the 5' *)
	begingn: int; (* begining of the segmented projection in the virtual cdna of the gene *)
	endingn: int; (* end of the segmented projection in the virtual cdna of the gene *)
	
 }

  let create chr gb ge st sc ep cov no bgn egn = { chrom=chr; gbeg=gb; gend=ge; str=st; score=sc; exproj=ep; coverage=cov; noinexproj=no; begingn=bgn; endingn=egn}
  
  let setcoverage cov s = {s with coverage=cov}
  let setscore sc s = {s with score=sc}

  let chrom s = s.chrom
  let gbeg s = s.gbeg
  let gend s = s.gend
  let str s = s.str
  let score s = s.score
  let exproj s = s.exproj
  let coverage s = s.coverage
  let noinexproj s = s.noinexproj
  let begingn s = s.begingn
  let endingn s = s.endingn
end



