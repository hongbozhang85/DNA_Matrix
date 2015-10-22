compiler = gfortran
flag = 

driver = driver.f90
driver1 = imagedriver.f90
driver2 = motscurve.f90
driver3 = strscurve.f90
objects = readdata.o type.o utils.o main.o generate_seq.o output.o cluster.o

help:
	@echo "genomics  : general analysis of gemone."
	@echo "          such as: cal S(l_s), generate 2-dim image, clustered image, "
	@echo "                   motif count histogram, curve S(l_m), shanon box
	@echo "image     : entropy curve from length=2 to 13, 2-dim image from length=2 to 13"
	@echo "            length=13 needs 768MB at least."
	@echo "motscurve : S(l_s,l_m,n_c)
	@echo "strscurve : pixel S from length=2 to a 13. include random."
	@echo "others...debugging..."

utils.o: type.o
readdata.o: type.o utils.o
generate_seq.o: type.o
main.o: type.o output.o
cluster.o: type.o
output.o: type.o generate_seq.o utils.o

%.o: %.f90
	$(compiler) $(flag) -c $*.f90

genomics: $(objects) $(driver)
	$(compiler) $(flag) $(objects) $(driver) isort.a -o $@
	cp genomics ../genomics 

image: $(objects) $(driver1)
	$(compiler) $(flag) $(objects) $(driver1) isort.a -o $@
	cp image ../image
	
motscurve: $(objects) $(driver2)
	$(compiler) $(flag) $(objects) $(driver2) isort.a -o $@
	cp motscurve ../motscurve
	
strscurve: $(objects) $(driver3)
	$(compiler) $(flag) $(objects) $(driver3) isort.a -o $@
	cp strscurve ../strscurve
	
clean: 
	rm -f *.o core *.mod
