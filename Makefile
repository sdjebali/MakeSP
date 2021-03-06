#######################
#    OCAML Programs   #
#######################

OCAMLC = ocamlc
OCAMLOPT = ocamlopt
OCAMLDEP = ocamldep



#######################
#     OCAML Flags     #
#######################

INCLUDES =					#all relevant -I options here
OCAMLFLAGS = $(INCLUDES)	#add other options for ocamlc here
OCAMLOPTFLAGS = $(INCLUDES)	#add other options for ocamlopt here



#######################
#    Sources files    #
#######################

MLFILES = common.ml config.ml collection.ml bio_objects.ml conversionES.ml input.ml compute_exproj.ml output.ml main.ml
CMXAFILES = $(CMAFILES:%.cma=%.cmxa) 
CMOFILES = $(MLFILES:%.ml=%.cmo) 
CMXFILES = $(MLFILES:%.ml=%.cmx)
CMIFILES = $(MLFILES:%.ml=%.cmi) $(MLIFILES:%.mli=%.cmi)
OBJFILES = $(CMIFILES:%.cmi=%.o)
BINFILE = makeSP




#######################
#        Rules        #
#######################

makeSP: $(CMXFILES)
	@echo "LNKOPT $(BINFILE)"
	@$(OCAMLOPT) $(CMXAFILES) $(CMXFILES) $(OCAMLOPTFLAGS) -o $(BINFILE)


# Common rules
.SUFFIXES: .ml .cmo .cmi .cmx

.ml.cmo:
	@echo "OCAMLC $<"
	@$(OCAMLC) $(OCAMLFLAGS) -c $<

.mli.cmi:
	@echo "OCAMLC $<"
	@$(OCAMLC) $(OCAMLFLAGS) -c $<

.ml.cmx:
	@echo "OCAMLOPT $<"
	@$(OCAMLOPT) $(OCAMLOPTFLAGS) -c $<

# Clean up
clean:
	@echo "Cleaning .cmo .cmx .o .cmi and binary"
	@rm -f $(BINFILE) $(CMOFILES) $(CMXFILES) $(CMIFILES) $(OBJFILES) 

# Dependencies
depend:
	@echo "Calculating dependencies"
	@$(OCAMLDEP) $(INCLUDES) $(MLFILES) $(MLIFILES) > $(DEPENDFILE)

# include $(DEPENDFILE)


