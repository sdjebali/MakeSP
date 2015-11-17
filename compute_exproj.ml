(* compute_exproj.ml *)

open Common
open Collection
open Bio_objects



(* cut a list of exons order by gene and then beg and end positions by gene *)
let cutintogenes sexons =
  let auxcut tex =
    let i = ref 1 and lpos = ref [0] and curgn = ref (Exon.gnid tex.(0)) and lgarr = Array.length tex in
      while (!i<lgarr) do
	if ((Exon.gnid tex.(!i)) <> !curgn) then
	  begin
	    lpos := !i::(!lpos);
	    curgn := Exon.gnid tex.(!i);
	  end;
	incr i;
      done;
      List.rev (!lpos) in 
    SegSeq.setsegment sexons (auxcut (SegSeq.tank sexons))



(* makes the assumption that the array included in sexofagene is ordered according to 
   the genomic begining and end positions.
 
   It would be nice this function could be generalised taking as parameter the functions 
   useful for the cutting , eg here Exon.gbeg and Exon.gend 
   like 
   let cutintoexproj sexofagene fbeg fend 
   instead of 
   let cutintoexproj sexofagene
*)
let cutintoexproj sexofagene =
  let auxcut tex =
    let i = ref 0 and lpos = ref [0] and curend = ref 0 and lgarr = Array.length tex in
      while (!i<lgarr) do
	curend := Exon.gend tex.(!i);
	
	while (!i<lgarr && ((Exon.gbeg tex.(!i)) <= !curend)) do
	  if ((Exon.gend tex.(!i)) > !curend) then
	    curend:= Exon.gend tex.(!i);
	  incr i;
	done;

	if(!i<lgarr) then
	  lpos := !i::(!lpos);

      done;
      List.rev (!lpos) in 
    SegSeq.setsegment sexofagene (auxcut (SegSeq.tank sexofagene))




(* Takes as input a list of arrays of exons of a gene, each array recording the exons of a 
   particular exon projection, and creates the corresponding list of exon projections 
*)
let exarrl2exprojl exarrl  =
  let cmlgth = ref 0 and exprojl = ref [] and i = ref 0 and n = List.length exarrl and fstarr = List.hd exarrl in
  let chrom = Exon.chrom fstarr.(0) and strand = Exon.str fstarr.(0) in

    while(!i<n) do
      (* the ith array of the list, corresponding to the ith exon projection *)
      let ai= List.nth exarrl !i in
	(* nb exon in the exon projection *)
      let ni = Array.length ai and exendi = Array.map Exon.gend ai in
	Array.sort Pervasives.compare exendi;
	(* not optimal since we have already computed the end (curend = max of the ends) in cutintoexproj! *)
	let gbeg = Exon.gbeg ai.(0) and gend = exendi.(ni-1) in
	  exprojl:= (ExonProj.create chrom gbeg gend strand ni ai (!i+1) (!cmlgth+1) (!cmlgth+gend-gbeg+1))::(!exprojl);
	  cmlgth:=!cmlgth+gend-gbeg+1;
	  incr i;
    done;
    List.rev (!exprojl)





(* this function is used by beglendl2intervals, this is here to remember which position
   was an exon begining position and which position was an exon end position in order
   that we get all sp, even of 1 nuclotide, and that none of them are overlapping by 1 base.
   intervlist is a simple list of the following type: [(a,b),(b,c),(c,d)...] but here we want 
   to adjust either the begining or the end or both of these intervals to take into account 
   the fact that these positions are begining or end of exons.

   The idea is to adjust each intervall of intervlist and also to add sp of size 1 nucleotide
   in case we have positions that are both in begexlist and in endexlist

   Note: it is not necessary to order the intervals as it will be done by the function beglendl2intervals
   However here we remove the redundancy that may be in the list of intervals 
*)
let find_all_sp_intervals begexlist endexlist intervlist =
  (* filter the elements of begexlist that are in endexlist *)
  let intervlist_adjusted = List.sort Pervasives.compare 
    (List.map (fun (b,e) -> 
		 if (b=e) then
		   (b,e)
		 else
		   begin
		     if (List.mem b endexlist) then 
		       begin
			 if (List.mem e begexlist) then
			   (* both b in endexlist and e in begexlist *)
			   begin
			     (b+1,e-1)   
			   end
			 else
			   (* only b in endexlist *)
			   begin
			     (b+1,e)  
			   end
		       end
		     else
		       begin
			 if (List.mem e begexlist) then   
			   (* only e in begexlist *)
			   begin
			     (b,e-1)    
			   end
			 else
			   (* neither b in endexlist nor e in begexlist *)
			   (b,e)
		       end
		   end
	      ) intervlist) in
    remove_redund Pervasives.compare intervlist_adjusted



(* Takes in a list of begining positions and a list of end positions that are supposed to be
   ORDERED, and computes the list of intervals, represented by these positions.
   minsz is the minimum accepted size/length for a segment/intervalle (eg 2 nt, which is the 
   minimum minimum size, see underneath)

   Important remark: begl and endl are ordered and without redundancy
*)
let beglendl2intervals begl endl minsz =
  (* merges the begining and end positions and order them.
     bE CAREFUL: here we loose track of which position was an exon begining position
     and which position was an exon end position.

     This has important consequences:
     1) The begining position of a segmented projection is equal to the end of the following 
     segmented projection when the segmented projection is not the first of the exon projection considered
     
     2) In such a case we do not know if it is the end position of the preceeding segmented projection that 
     has to be decreased by 1 or the begining position of the current segmented projection that has to be 
     increased by 1, since this depends on the fact that the position was a begining position or an end 
     position of a particular exon!!! (and what if it is both??)
 
     3) The minimum length tolerated for a segmented projection is 2, which in fact means 1 for all segmented 
     projections except the first one of the exon projection (--> minor)

     We need to compute the right begining positions for non initial segmented projections
     *************************************************************************************
     In fact we need to take as input of our interval making function the list of begining and the list of end positions
     and maybe make an iterative function instead of a recursive one...

     Or if not to compute the intervals as we do it now and then to modify the list in the following way:
     ****************************************************************************************************
     If the begining position of such an interval is a begining exon position and not an end exon position:
           * keep it like it is, if not increase it by 1
     ElseIf the end position of such an interval is an end exon position and not a begining exon position: 
           * keep it like this, if not decrease if by 1
     Else create and add an interval of 1 nucleotide corresponding to THIS position and also:
     1) decrease the interval end position which is equal to this position
     2) increase the interval begining position which is equal to this position

     The intervals that are output by this function are ordered according to begining position
     *****************************************************************************************
  *)
  let begendl = List.merge Pervasives.compare begl endl in
    List.sort (fun i1 i2 -> Pervasives.compare (fst i1) (fst i2)) (List.filter (fun (a,b) -> b-a+1>=minsz) (find_all_sp_intervals begl endl (Common.intervals begendl)))
      





(* Note: this is an initial creation for the segmented projection, it assigns the coverage and the score to 0.
   The correct coverage and score will be computed later on. 
   It will be good to know the beg and end positions in gn of each SP
*)
let interv2segp chrom strand exproj i (a,b)=
  (* i+1 is the number of the segmented projection in the exon projection *)
  SegProj.create chrom a b strand 0.0 exproj 0 (i+1) 0 0


(* Be careful the ends of the exon are NOT ordered *)
let exp2segparr minsz exproj =
  (* exons, chromosome and strand of the exon projection *)
  let exarr = ExonProj.exarr exproj and chrom = ExonProj.chrom exproj and strand = ExonProj.str exproj in
    (* begining and end positions of the exons composing the exon projection, these two lists are ORDERED and without redundancy *)
  let begl = remove_redund Pervasives.compare (List.sort Pervasives.compare (Array.to_list (Array.map Exon.gbeg exarr))) and endl = remove_redund Pervasives.compare (List.sort Pervasives.compare  (Array.to_list (Array.map Exon.gend exarr))) in
    (* for each interval make a segmented projection *)
    Array.mapi (interv2segp chrom strand exproj) (Array.of_list (beglendl2intervals begl endl minsz))
  


(* To compute the coverage of each segmented projection we need to know to how many exons a given nucleotide
   of the segmented projection (let us say the 3rd one) belongs to *)
let poswithin2 pos (b,e,_) =
  b<=pos && pos<=e


(* simplest scoring function for computing the score of a given segmented projection
   given the list of exons (or rather the list of triplets (beg,end,score) of the exons)
   it is composed of: ie adds the scores of all its components *) 
let rec addingfunction begendsclist =
  match begendsclist with
    |[] -> 0.0
    |(b,e,sc)::q -> sc+.(addingfunction q)


(* nth is the number of the nucleotide of the segmented projection (counting from the begining 
   and departing from 1) that will be tested for the coverage, usually 3 but has to do with minsz above.
   scorefunction is a function that compute the score of a segmented projection given the list of triplets
   (begpos, endpos, score) of each exon this segmented proj is composed of.
*)
let segp2segpwithcov scorefunction nth segp =
  let exlist = Array.to_list (ExonProj.exarr (SegProj.exproj segp)) in
  let begendsclistofexproj = List.rev (List.rev_map (fun ex -> ((Exon.gbeg ex), (Exon.gend ex), (Exon.score ex))) exlist) in
  let begendsclistofsproj = List.filter (poswithin2 ((SegProj.gbeg segp)+nth-1)) begendsclistofexproj in
   (* the coverage is the number of transcripts a segmented projection is included in *)
  let cov = List.length begendsclistofsproj in
  let sc = scorefunction begendsclistofsproj in  
    SegProj.setscore sc (SegProj.setcoverage cov segp)
  

