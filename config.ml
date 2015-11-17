(* config.ml *)

(*********************************************************)
(* Structure for the parameters of the program makeSP    *)
(*********************************************************)
type 'a context_t =
	{ 
	  mutable file:string;			(* input file *)
	  mutable outformat: string;            (* format of output file, either txt or gff with all the info, 
						   or bed with only the segmented projections and the coverage *)
	  mutable verbose: bool;                (* if we want the output to be displayed on stdout *)
	} 



(********************************************************)
(* overlap context, these are the default parameters    *)
(********************************************************)
let context = 
  {	
    file = "";
    outformat="gff";
    verbose=false;
  };;


let usage = 
" 
                     *********************************************          
                     *   makeSP - version v1.0 (February 2009)   *
                     *             Sarah Djebali                 *
                     *********************************************

Usage: "^(Filename.basename (Sys.argv.(0)))^" file [options] 

Takes a set of exons associated to transcript and gene names and does two things:
1) project all the exons of each gene on the genomic sequences
2) identifies for each such projected exon of a gene, all the segmented projections
it is composed of, i.e. all the maximal segments having the same transcript coverage.

** file must be provided in version 2 gff format. Moreover if needs to fulfill the following conditions:
   1) the 4 first subfields of the 9th field must represent the transcript_id and the gene_id of the exon 
   2) the transcript and the gene ids must be surrounded by double quotes and be followed by a semicolon
   here is an example:
   chr7    VEGA_Known      exon    115444498       115444739       .       +       .       transcript_id \"AC073130.2-001\"; gene_id \"TES\"; gene_alias \"AC073130.2\"; exon_id \"AC073130.2-001-exon1\";

** [options] can be:
   -f outformat:  outformat is the format of output file the user wants to get. It can either be txt, gff or bed.
                  Note: the bed format is more compact but only contains information about the segmented projections and the coverage.
                  Note: the name of the output file will be [file]_segproj.[outformat]
                  -> default is gff

   -v:            displays the result on the standard output.
                  -> default is false

** Please report any bug to sarahqd@gmail.com
"

(*
   -st strmode:    defines whether and how strand is taken into account when computing overlap between two features:
                  * strmode=0: strand is not taken into account, only positions are
                  * strmode=1: only features on the same strand can be included in the same projection (unstranded features can match anything)  *** NOT implemented yet
                  -> default is 0

*)
(***********************************************************************)
(* Read the arguments from the command line and updates the context    *)
(***********************************************************************)
let read_commandline () =
  let _u = try 
      begin
	context.file <- Sys.argv.(1);
      end
    with
      | Invalid_argument s -> Common.print_error usage
  in
  
  (* we start reading the arguments from the 2nd one since the first is compulsory and is the input file *)
  let argnum = ref 1 and ok = ref true in

  (* This function returns the next argument. mustbe says whether this argument has to exist.
     In general mustbe is true since the arguments go by pairs *)
  let getarg mustbe =
    incr argnum; 
    if !argnum < (Array.length Sys.argv) then 
      Sys.argv.(!argnum)
    else 
      if mustbe then 
	raise Not_found 
    else
      ""
  in
    (* Actually reading each of the arguments *)
    try 
      while (!ok) do
	match (getarg false) with
	  | "-f"  -> context.outformat <- getarg true
 (*	  | "-st"  -> context.strmode <- int_of_string (getarg true) *)
	  | "-v"  -> context.verbose <- true
	  | "-h"
	  | "--help" -> Printf.printf "%s\n" usage; exit 1; 
	  | ""	-> ok := false
	  | s	-> failwith s
      done;
      Common.print_log "# Command line has been read\n";
    with
      | Not_found -> Common.print_error ("Missing parameter at the end of the command line\n"^usage^"\n");
      | Failure "int_of_string" -> Common.print_error ("Syntax error (incorrect integer value)\n"^usage^"\n");
      | Failure s -> Common.print_error ("Syntax error ("^s^")\n"^usage^"\n");;


