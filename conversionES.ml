(* conversionES.ml *)

open Common
open Bio_objects



let strand_of_string s =
  match s with
    | "++"
    | "+" -> Forward
    | "+-"
    | "-" -> Reverse
    |  _  -> print_error "The strand should be provided as the 7th field of your gtf file\n" ; raise Not_found;;


let strand_to_string = function
  | Forward -> "+"
  | Reverse -> "-";;


let strand_to_string2 = function
  | Forward -> "p"
  | Reverse -> "m";;







