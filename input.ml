(* input.ml *)
(* reads input files and store the data in appropriate structures *)


open Common
open Bio_objects
open ConversionES



(*
  44regions_CHR_coord_exons_Vegasure.gtf
  chr7    VEGA_Known      exon    115444498       115444739       .       +       .       transcript_id "AC073130.2-001"; gene_id "TES"; gene_alias "AC073130.2"; exon_id "AC073130.2-001-exon1";
*)


(* record_of_line_gff takes as input a line of a gff file (lline) and outputs an exon. 
   Note : to extend this function to be able to read other gff files than including exons
   we could make use of tab.(6) 
*)
let record_of_line_gff lline =
  let tab = Array.of_list lline in
  let n = Array.length tab in
    if (n<=8) then
       begin
	failwith "record_of_line_gff: Syntax error in your gff file: it should contain at least 9 fields separated by tabs"
       end
    else
      begin
	let attl = split ' ' (Common.clean_beg_string (Common.clean_end_string tab.(8))) in
	let natt = List.length attl in
	  if((natt mod 2)!=0) then
	    begin
	      failwith "record_of_line_gff: Syntax error in your gff file: it should contain at least 9 fields separated by tabs, and the 9th field should be divided into a set of attibutes of the form [key space value] separated by spaces. Also there should not be any space at the end"
	    end
	  else
	    begin
	      if (natt <= 3) then
		failwith "record_of_line_gff: Syntax error in your gff file: it should contain at least 9 fields separated by tabs, and the 9th field should be divided into a set of attibutes of the form [key space value] separated by spaces, the two first attibutes being the transcript_id and the gene_id. Be careful transcript and gene ids must be surrounded by double quotes and followed by a semicolon."
	      else
		begin
		  let name_tr_tmp = List.nth attl 1 and name_gn_tmp = List.nth attl 3 in
		  let name_tr = try (String.sub name_tr_tmp 1 ((String.length name_tr_tmp) - 3)) with
		      Invalid_argument _ -> failwith "record_of_line_gff: Syntax error in your gff file: it should contain at least 9 fields separated by tabs, and the 9th field should be divided into a set of attibutes of the form [key space value] separated by spaces, the two first attibutes being the transcript_id and the gene_id. Be careful transcript and gene ids must be surrounded by double quotes and followed by a semicolon."; ""
		  and name_gn = try (String.sub name_gn_tmp 1 ((String.length name_gn_tmp) - 3)) with
		      Invalid_argument _ -> failwith "record_of_line_gff: Syntax error in your gff file: it should contain at least 9 fields separated by tabs, and the 9th field should be divided into a set of attibutes of the form [key space value] separated by spaces, the two first attibutes being the transcript_id and the gene_id. Be careful transcript and gene ids must be surrounded by double quotes and followed by a semicolon."; ""
		  in
   		    Exon.create 
		      tab.(0)
		      (int_of_string tab.(3))
		      (int_of_string tab.(4))
		      (strand_of_string tab.(6))
		      (try (float_of_string tab.(5)) with Failure s -> 0.0)
		      tab.(2)
		      name_tr
		      name_gn
		end
	    end
      end


(* This function reads the lines of a (version 2) gff exons file (input channel), 
   and outputs the corresponding list of exons. 
   Note: still too slow on big files.
*)
let make_exon_list inchan =   
  let l = ref [] and stop = ref false in
    while (not !stop) do
      (
	try
	  let currline = (split '\t' (input_line inchan)) in
	  let currexon= record_of_line_gff currline in
	    l:=(currexon)::(!l)
	with
	  | End_of_file -> stop:=true
      )
    done;
    !l;;
    

