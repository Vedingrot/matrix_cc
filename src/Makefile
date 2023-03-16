CC          = g++ $(STANDART) $(CFLAGS)
STANDART		= --std=c++17
CDEBUG      = -g
CFLAGS      = $(CDEGUB) -Wall -Werror -Wextra -Weffc++
HEADER      = s21_matrix_oop.h
SOURCES     = s21_matrix_oop.cc
OBJECTS     = $(SOURCES:.cc=.o)
LIB_MATRIX  = libs21_matrix_oop.a
TEST_SOURCE 				= test.cc
TEST_FLAGS  			  = -lgtest
ifeq ($(shell uname), Linux)
	TEST_FLAGS += -lpthread
endif
TEST_OUT    = test.out
RM          = rm -rf
GCOV_COMPILE_FLAGS  = -fprofile-arcs -ftest-coverage
GCOV_OUT		= gcov_report.out
GCOV_RESULT	= *.gcda *.gcno *.gcov
GCOV 			  = gcov
GCOV_FLAGS  = -kjr

.PHONY: all
all: $(LIB_MATRIX)

ARCHIVATOR = ar rcs
$(LIB_MATRIX): $(OBJECTS)
	$(ARCHIVATOR) $(LIB_MATRIX) $^

.PHONY: clean
clean:
	$(RM) $(TEST_OUT) $(LIB_MATRIX) $(OBJECTS) $(GCOV_OUT) $(GCOV_RESULT)\
		$(LCOV_REPORT_DIR) $(COVERAGE_INFO)

$(TEST_OUT): $(TEST_SOURCE) $(LIB_MATRIX)
	$(CC) $(CFLAGS) $^ $(TEST_FLAGS) -o $@

test: $(TEST_OUT)
	./$(TEST_OUT)

$(GCOV_OUT): $(TEST_SOURCE) $(SOURCES)
	$(RM) $(GCOV_RESULT)
	$(CC) $^ $(TEST_FLAGS) $(GCOV_COMPILE_FLAGS) -o $@

gcov_report: $(GCOV_OUT)
	./$(GCOV_OUT)
	$(GCOV) $(GCOV_FLAGS) $(SOURCES)

ifeq ($(shell uname), Darwin)
BREW_PATH += ~/.brew/bin/
endif
LCOV          = $(BREW_PATH)lcov
LCOV_FLAGS    = --no-external -c -d .
COVERAGE_INFO = coverage.info

$(COVERAGE_INFO): $(GCOV_OUT)
	./$(GCOV_OUT)
	$(LCOV) $(LCOV_FLAGS) -o $(COVERAGE_INFO)

LCOV_REPORT_DIR = lcov_report_dir
GENHTML         = $(BREW_PATH)genhtml

$(LCOV_REPORT_DIR) : $(COVERAGE_INFO)
	$(GENHTML) $^ -o $(LCOV_REPORT_DIR)

lcov_report: $(LCOV_REPORT_DIR)

$(OBJECTS): $(SOURCES)
	$(CC) -c $(SOURCES)

$(LIB_MATRIX): $(HEADER) $(SOURCES)

.IGNORE: linter cppcheck style check

.PHONY: style_check cppcheck linter leaks
style_check: cppcheck linter

SUPPRESSING = --suppress=unusedFunction --suppress=missingIncludeSystem \
			  --suppress=unmatchedSuppression --suppress=missingInclude
cppcheck:
	cppcheck --enable=all $(SUPPRESSING) $(STANDART) $(SOURCES)

linter:
	cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	python3 ../materials/linters/cpplint.py --extensions=cc *.h *.cc
	rm -f CPPLINT.cfg

leaks: $(TEST_OUT)
	leaks --atExit -- ./$^
