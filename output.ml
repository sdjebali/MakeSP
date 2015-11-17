(* output.ml *)


open Bio_objects
open ConversionES


(* o is the output channel corresponding to the file where we need to write the segmented projection
   Examples:
   - gff format
   chr1    makeSP  segproj 8       8       0.000000        -       .       exproj: chr1_8_216_m gene: A noep_ingn: 1 begep_ingn: 1 endep_ingn: 209 spcoverage: 1 nosp_inep: 1

   - txt format
   #SPname    EPname    GNname    NoEPinGN    BegEPinGN    EndEPinGN    SPcoverage    NoSPinEP  SPscore
   SegProj_chr21_14568284_14568349_p       ExonProj_chr21_14568284_14568349_p      ABCC13  1       1       66      1       1

   - bed format
   chr22   14430401        14430427        A       1       +


   The difference with the first version is that we have added the score of the segmented projection at the end
*)
let print o format segp =
  let exp = SegProj.exproj segp in
  let exarr = ExonProj.exarr exp in
    match format with

      | "gff" -> Printf.fprintf o "%s\tmakeSP\tsegproj\t%i\t%i\t%f\t%s\t.\texproj: %s_%i_%i_%s gene: %s noep_ingn: %i begep_ingn: %i endep_ingn: %i spcoverage: %i nosp_inep: %i\n" (SegProj.chrom segp) (SegProj.gbeg segp) (SegProj.gend segp) (SegProj.score segp) (strand_to_string (SegProj.str segp)) (ExonProj.chrom exp) (ExonProj.gbeg exp) (ExonProj.gend exp) (strand_to_string2 (ExonProj.str exp)) (Exon.gnid exarr.(0)) (ExonProj.noingn exp) (ExonProj.begingn exp) (ExonProj.endingn exp) (SegProj.coverage segp) (SegProj.noinexproj segp)

      | "txt" -> Printf.fprintf o "SegProj_%s_%i_%i_%s\tExonProj_%s_%i_%i_%s\t%s\t%i\t%i\t%i\t%i\t%i\t%f\n" (SegProj.chrom segp) (SegProj.gbeg segp) (SegProj.gend segp) (strand_to_string2 (SegProj.str segp)) (ExonProj.chrom exp) (ExonProj.gbeg exp) (ExonProj.gend exp) (strand_to_string2 (ExonProj.str exp)) (Exon.gnid exarr.(0)) (ExonProj.noingn exp) (ExonProj.begingn exp) (ExonProj.endingn exp) (SegProj.coverage segp) (SegProj.noinexproj segp) (SegProj.score segp) 

      | "bed" -> Printf.fprintf o "%s\t%i\t%i\t%s\t%i\t%s\n" (SegProj.chrom segp) ((SegProj.gbeg segp)-1) (SegProj.gend segp) (Exon.gnid exarr.(0)) (SegProj.coverage segp) (strand_to_string (ExonProj.str exp)) (* theoretically the strand of a segmented projection can be undetermined, in case it is made of exons on different strands. This is not taken care of still. *)

      | _ -> raise (Invalid_argument ("The output format must be gff, txt or bed"));;
    

 
