# GPU backend for SDPB (based on SDPviBe by Antoine Tilloy, the prototype; see
# README.md and docs/GPU_backend_credits.md): standalone build of the device
# library and its checks, beside the CPU build (waf integration: later stage).
#
#   module load cuda/12.9.0   # for example; any CUDA >= 12.2
#   make GMP_DIR=<prefix of GMP> [ARCH="90 89"]
#
NVCC    ?= nvcc
ARCH    ?= 90
GMP_DIR ?= /usr
MPFR_DIR ?= $(GMP_DIR)
BUILD   ?= build

GENCODE := $(foreach a,$(ARCH),-gencode arch=compute_$(a),code=sm_$(a))
NVFLAGS := -std=c++17 -O3 $(GENCODE) -Xcompiler -Wall -Xcompiler -fopenmp -I$(GMP_DIR)/include
LDLIBS  := -lcublasLt -L$(GMP_DIR)/lib -Xlinker -rpath -Xlinker $(GMP_DIR)/lib -lgmp -lopenblas -lgomp

LIB_SRC := rns/Rns_Basis.cxx rns/Gram_Accumulator.cu rns/Mpf_Gram.cu rns/Rns_Gemm.cu linalg/Block_Factor.cu solver/Mpf_Ops.cu solver/Gpu_Solver.cu
LIB_OBJ := $(patsubst %,$(BUILD)/%.o,$(LIB_SRC))

all: $(BUILD)/gram_check $(BUILD)/mpf_gram_check $(BUILD)/gemm_check $(BUILD)/chol_tile_check $(BUILD)/mpf_gemm_check $(BUILD)/block_factor_check $(BUILD)/real_block_check $(BUILD)/q_check $(BUILD)/residual_check $(BUILD)/mpf_ops_check $(BUILD)/factor_profile $(BUILD)/gemv_compare

$(BUILD)/%.cxx.o: %.cxx $(wildcard rns/*.hxx rns/*.cuh linalg/*.hxx linalg/*.cuh solver/*.hxx)
	@mkdir -p $(dir $@)
	$(NVCC) $(NVFLAGS) -x cu -c $< -o $@

$(BUILD)/%.cu.o: %.cu $(wildcard rns/*.hxx rns/*.cuh linalg/*.hxx linalg/*.cuh solver/*.hxx)
	@mkdir -p $(dir $@)
	$(NVCC) $(NVFLAGS) -c $< -o $@

$(BUILD)/gram_check: $(BUILD)/test/gram_check.cu.o $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/mpf_gram_check: $(BUILD)/test/mpf_gram_check.cu.o $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/gemm_check: $(BUILD)/test/gemm_check.cu.o $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/mpf_gemm_check: $(BUILD)/test/mpf_gemm_check.cu.o $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/block_factor_check: $(BUILD)/test/block_factor_check.cu.o $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/real_block_check: $(BUILD)/test/real_block_check.cu.o $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/q_check: test/q_check.cu $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) -Xcompiler -fopenmp $^ -o $@ $(LDLIBS)

$(BUILD)/mpf_ops_check: $(BUILD)/test/mpf_ops_check.cu.o $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/gemv_compare: $(BUILD)/test/gemv_compare.cu.o $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/residual_check: $(BUILD)/test/residual_check.cu.o $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/factor_profile: $(BUILD)/test/factor_profile.cu.o $(LIB_OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/ref_solve: test/ref_solve.c
	@mkdir -p $(BUILD)
	$(CC) -O3 -fopenmp -I$(GMP_DIR)/include -I$(MPFR_DIR)/include $< -o $@ -L$(GMP_DIR)/lib -L$(MPFR_DIR)/lib -Wl,-rpath,$(GMP_DIR)/lib -Wl,-rpath,$(MPFR_DIR)/lib -lmpfr -lgmp -lm

$(BUILD)/compare_references: test/compare_references.c
	@mkdir -p $(BUILD)
	$(CC) -O2 -I$(GMP_DIR)/include -I$(MPFR_DIR)/include $< -o $@ -L$(GMP_DIR)/lib -L$(MPFR_DIR)/lib -Wl,-rpath,$(GMP_DIR)/lib -Wl,-rpath,$(MPFR_DIR)/lib -lmpfr -lgmp -lm

mpfr-checks: $(BUILD)/ref_solve $(BUILD)/compare_references

$(BUILD)/chol_tile_check: test/chol_tile_check.cu linalg/basecase.cuh
	@mkdir -p $(BUILD)
	$(NVCC) $(NVFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -rf $(BUILD)

.PHONY: all clean mpfr-checks
