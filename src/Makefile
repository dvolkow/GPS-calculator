#------------------------------------------
dirl := $(PATH) #--собственно папка где лежат модули
swdk := $(addsuffix /*.o,$(dirl))   #--добавление пути к файлам 
wwdk := $(addsuffix /*.f03,$(dirl))
comp := gfortran  #--компилятор
opt  := -o #--опция компилятора для исполняемой программы
op   := -O2 -ffree-line-length-none      #--параметры оптимизации
obh  := -c        #--опция компилятора для объектных модулей
#-правило по умолчанию
all: $(patsubst %.f03,%.o,$(wildcard *.f03)) $(patsubst %.f03,%.o,$(notdir $(wildcard $(wwdk))))
	$(comp) $(op) $^ $(opt) ../gpsCalc
	rm *.o
%.o : %.f03
	$(comp) $(obh) $(op) $<
