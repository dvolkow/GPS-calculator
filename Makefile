#VPATH := pwd #--объявление директории где лежат модули
#------------------------------------------
dirl := pwd   #--собственно папка где лежат модули
swdk := $(addsuffix /*.o,$(dirl))   #--добавление пути к файлам 
wwdk := $(addsuffix /*.f03,$(dirl))
comp := gfortran  #--компилятор
opt  := -o        #--опция компилятора для исполняемой программы
op   := -Ofast       #--параметры отпимизации
obh  := -c        #--опция компилятора для объектных модулей
#some_file :=  #--ручное перечисление нужных файлов для второго правила
#-правило по умолчанию, когда нет лишнего
all: $(patsubst %.f03,%.o,$(wildcard *.f03)) $(patsubst %.f03,%.o,$(notdir $(wildcard $(wwdk))))
	$(comp) $(op) $^ $(opt) programma
%.o : %.f03
	$(comp) $(obh) $(op) $<
#-удаление объектных модулей
.PHONY : clean
clean: 
	rm *.o
#-запуск программы
exe:
	./programma