(* main.ml: compute the segmented projections of a set of gene with multiple transcripts
   provided in a gtf file 
   VERY important note: the projection is done within each gene
*)

open Common
open Config
open Collection
open Bio_objects
open ConversionES
open Input
open Compute_exproj
open Output

(* starts from a gtf file like: 44regions_CHR_coord_exons_Vegasure.gtf
   compute the segmented projections, their coverage and their score
   taking into account the score of the exons constituting the segproj.

   Be careful: there are two parameters that have not been taken out here: 
   - the integer 3 used in the exp2segprr function
     This is the minimum length of a segmented projection for it to be output
     by the program.
   - the integer 2 used in the segp2segpwithcov function and which is the nucleotide
     that will be tested to compute the transcript coverage (ref 1). It depends on the 
     preceding parameter.
 *)




let computeSP () =
  read_commandline ();

  let earr = Array.of_list (Input.make_exon_list (open_in context.file)) in
  let _u = Common.print_log "# I have read the gtf file\n" in

  (* Orders the exons by gene and then by begining and end positions *)
  let _u = Array.sort Exon.compare earr in
  let _u = Common.print_log "# I have sorted the exons\n" in

  (* Cuts the set of exons according to the gene they belong to *)
  let s = cutintogenes (SegSeq.make2 earr) in
  let _u = Common.print_log "# I have separated the exons according to the genes\n" in
  
  (* Produces a list of exon array, one exon array for each gene *)
  let ltex= SegSeq.elements s in
  let ls = List.rev (List.rev_map SegSeq.make2 ltex) in
  let _u = Common.print_log "# I have produced one exon array per gene\n" in
    
  (* For each gene cuts the exon set into sets that will each compose an exon projection *)
  let ls2 = List.rev (List.rev_map cutintoexproj ls) in
  let lla = List.rev (List.rev_map SegSeq.elements ls2) in
  let llexproj = List.rev (List.rev_map exarrl2exprojl lla) in
  let _u = Common.print_log "# I have produced the projected exons for each gene\n" in
    
  (* Makes a set of segmented projections from each exon projection - Makes a segmentation fault for big files 
     and this is due to one of the List.map calls -> transform to List.rev (List.rev_map 
  *)
  let lltsegs = List.rev (List.rev_map (fun lexp -> List.rev (List.rev_map (fun exp -> exp2segparr 1 exp) lexp)) llexproj) in
  let _u = Common.print_log "# I have produced the segmented projections of each projected exon\n" in
    
  (* Computes the coverage and the score of each segmented projection, this function could be made a an argument of the command line 
     I also had to transform List.map calls to List.rev (List.rev_map calls because of seg fault probably due to a stack overflow
 *)
  let lltsegswithcov = List.rev (List.rev_map (fun ltsegs -> List.rev (List.rev_map (fun tsegs -> Array.map (segp2segpwithcov addingfunction 1) tsegs) ltsegs)) lltsegs) in
  let _u = Common.print_log "# I have computed the coverage and score of each segmented projection\n" in
  let _u = Common.print_log "# I have done my job, I will now write the output file\n" in
  
    (* Writes the output in standard output or in context.outfile *)
    if(context.verbose) then
      begin
	List.iter (fun ltseg -> List.iter (fun tseg -> Array.iter (Output.print stdout context.outformat) tseg) ltseg) lltsegswithcov
      end
    else
      begin
	let o = open_out ((context.file)^("_segproj.")^(context.outformat)) in
	let _u = List.iter (fun ltseg -> List.iter (fun tseg -> Array.iter (Output.print o context.outformat) tseg) ltseg) lltsegswithcov in
	let _u = flush o in
	  ()
      end;;

computeSP ();;

