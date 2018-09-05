# Author: corbett@caltech.edu

from quintuple.computer import (Program)


class Programs(object):
    """Some useful programs collected in one place for running on the quantum computer class"""
    program_blue_state = Program("""h q[1];
            t q[1];
            h q[1];
            t q[1];
            h q[1];
            t q[1];
            s q[1];
            h q[1];
            t q[1];
            h q[1];
            t q[1];
            s q[1];
            h q[1];
            bloch q[1];""")

    program_test_XYZMeasureIdSdagTdag = Program("""sdg q[0];
            x q[1];
            x q[2];
            id q[3];
            z q[4];
            tdg q[0];
            y q[4];
            measure q[0];
            measure q[1];
            measure q[2];
            measure q[3];
            measure q[4];""")

    program_test_cnot = Program("""x q[1];
            cx q[1], q[2];""")

    program_test_many = Program("""sdg q[0];
            x q[1];
            x q[2];
            id q[3];
            z q[4];
            tdg q[0];
            cx q[1], q[2];
            y q[4];
            measure q[0];
            measure q[1];
            measure q[2];
            measure q[3];
            measure q[4];""")
    # IBM Tutorial Section III, Page 4
    program_zz = Program("""h q[1];
        cx q[1], q[2];
        measure q[1];
        measure q[2];""")  # "00",0.5; "11",0.5 # <zz> = 2

    program_zw = Program("""h q[1];
        cx q[1], q[2];
        s q[2];
        h q[2];
        t q[2];
        h q[2];
        measure q[1];
        measure q[2]""")  # "00",0.426777; "01",0.073223; "10",0.073223; "11",0.426777 # <zw> = 1/sqrt(2)

    program_zv = Program("""h q[1];
        cx q[1], q[2];
        s q[2];
        h q[2];
        tdg q[2];
        h q[2];
        measure q[1];
        measure q[2];""")  # "00",0.426777; "01",0.073223; "10",0.073223; "11",0.426777 # <zv> = 1/sqrt(2)

    program_xw = Program("""h q[1];
        cx q[1], q[2];
        h q[1];
        s q[2];
        h q[2];
        t q[2];
        h q[2];
        measure q[1];
        measure q[2];""")  # "00",0.426777; "01",0.073223; "10",0.073223; "11",0.426777 # <xw> =

    program_xv = Program("""h q[1];
        cx q[1], q[2];
        h q[1];
        s q[2];
        h q[2];
        tdg q[2];
        h q[2];
        measure q[1];
        measure q[2];""")  # "00",0.073223; "01",0.426777; "10",0.426777; "11",0.073223; # <xv> =

    # Currently not used, but creats a superposition of 00 and 01
    program_00_01_super = Program("""sdg q[1];
        t q[1];
        t q[1];
        s q[1];
        h q[1];
        h q[0];
        h q[1];
        h q[0];
        h q[1];
        cx q[0], q[1];
        measure q[0];
        measure q[1];""")
    # IBM Tutorial Section III, Page 5
    program_ghz = Program("""h q[0];
        h q[1];
        x q[2];
        cx q[1], q[2];
        cx q[0], q[2];
        h q[0];
        h q[1];
        h q[2];
        measure q[0];
        measure q[1];
        measure q[2];""", result_probability=[0.5, 0, 0, 0, 0, 0, 0, 0.5])  # "000":0.5; "111":0.5

    program_ghz_measure_yyx = Program("""h q[0];
        h q[1];
        x q[2];
        cx q[1], q[2];
        cx q[0], q[2];
        h q[0];
        h q[1];
        h q[2];
        sdg q[0];
        sdg q[1];
        h q[2];
        h q[0];
        h q[1];
        measure q[2];
        measure q[0];
        measure q[1];""", result_probability=[0.25, 0, 0, 0.25, 0, 0.25, 0.25,
                                              0])  # "000":0.25; "011": 0.25; "101": 0.25; "110":0.25

    program_ghz_measure_yxy = Program("""h q[0];
        h q[1];
        x q[2];
        cx q[1], q[2];
        cx q[0], q[2];
        h q[0];
        h q[1];
        h q[2];
        sdg q[0];
        h q[1];
        sdg q[2];
        h q[0];
        measure q[1];
        h q[2];
        measure q[0];
        measure q[2];""", result_probability=[0.25, 0, 0, 0.25, 0, 0.25, 0.25,
                                              0])  # "000":0.25; "011": 0.25; "101": 0.25; "110":0.25

    program_ghz_measure_xyy = Program("""h q[0];
        h q[1];
        x q[2];
        cx q[1], q[2];
        cx q[0], q[2];
        h q[0];
        h q[1];
        h q[2];
        h q[0];
        sdg q[1];
        sdg q[2];
        measure q[0];
        h q[1];
        h q[2];
        measure q[1];
        measure q[2];""", result_probability=[0.25, 0, 0, 0.25, 0, 0.25, 0.25,
                                              0])  # "000":0.25; "011": 0.25; "101": 0.25; "110":0.25

    program_ghz_measure_xxx = Program("""h q[0];
        h q[1];
        x q[2];
        cx q[1], q[2];
        cx q[0], q[2];
        h q[0];
        h q[1];
        h q[2];
        h q[0];
        h q[1];
        h q[2];
        measure q[0];
        measure q[1];
        measure q[2];""", result_probability=[0, 0.25, 0.25, 0, 0.25, 0, 0,
                                              0.25])  # "001":0.25; "010": 0.25; "100": 0.25; "111":0.25

    # IBM Tutorial Section IV, Page 1
    program_reverse_cnot = Program("""x q[2];
        h q[1];
        h q[2];
        cx q[1], q[2];
        h q[1];
        h q[2];
        measure q[1];
        measure q[2];""", result_probability=(0.0, 0.0, 0.0, 1.0))  # "11": 1.0

    program_swap = Program("""x q[2];
        cx q[1], q[2];
        h q[1];
        h q[2];
        cx q[1], q[2];
        h q[1];
        h q[2];
        cx q[1], q[2];
        measure q[1];
        measure q[2];""", result_probability=(0.0, 0.0, 1.0, 0.0))  # "10": 1.0

    program_swap_q0_q1 = Program("""h q[0];
        cx q[0], q[2];
        h q[0];
        h q[2];
        cx q[0], q[2];
        h q[0];
        h q[2];
        cx q[0], q[2];
        cx q[1], q[2];
        h q[1];
        h q[2];
        cx q[1], q[2];
        h q[1];
        h q[2];
        cx q[1], q[2];
        cx q[0], q[2];
        h q[0];
        h q[2];
        cx q[0], q[2];
        h q[0];
        h q[2];
        cx q[0], q[2];
        bloch q[0];
        bloch q[1];
        bloch q[2];""", bloch_vals=(
        (0, 0, 1), (1, 0, 0), (0, 0, 1), None, None))  # Bloch q0: (0,0,1); #q1: (1,0,0) q2: (0,0,1)
    program_controlled_hadamard = Program("""h q[1];
        s q[1];
        h q[2];
        sdg q[2];
        cx q[1], q[2];
        h q[2];
        t q[2];
        cx q[1], q[2];
        t q[2];
        h q[2];
        s q[2];
        x q[2];
        measure q[1];
        measure q[2];""", result_probability=[0.5, 0.0, 0.25, 0.25])  # "00": 0.5; "10": 0.25; "11":0.25
    program_approximate_sqrtT = Program("""h q[0];
        h q[1];
        h q[2];
        h q[3];
        h q[4];
        bloch q[0];
        h q[1];
        t q[2];
        s q[3];
        z q[4];
        t q[1];
        bloch q[2];
        bloch q[3];
        bloch q[4];
        h q[1];
        t q[1];
        h q[1];
        t q[1];
        s q[1];
        h q[1];
        t q[1];
        h q[1];
        t q[1];
        s q[1];
        h q[1];
        t q[1];
        h q[1];
        t q[1];
        h q[1];
        bloch q[1];""", bloch_vals=((1, 0, 0), (0.927, 0.375, 0.021), (0.707, 0.707, 0.000), (0.000, 1.000, 0.000), (
        -1.000, 0.000,
        0.000)))  # Bloch coords q0: (1.000, 0.000, 0.000) q1: (0.927, 0.375, 0.021) q2: (0.707, 0.707, 0.000) q3: (0.000, 1.000, 0.000) q4: (-1.000, 0.000, 0.000) # checks out when we manually get_bloch
    program_toffoli_state = Program("""h q[0];
        h q[1];
        h q[2];
        cx q[1], q[2];
        tdg q[2];
        cx q[0], q[2];
        t q[2];
        cx q[1], q[2];
        tdg q[2];
        cx q[0], q[2];
        t q[1];
        t q[2];
        cx q[1], q[2];
        h q[1];
        h q[2];
        cx q[1], q[2];
        h q[1];
        h q[2];
        cx q[1], q[2];
        cx q[0], q[2];
        t q[0];
        h q[1];
        tdg q[2];
        cx q[0], q[2];
        measure q[0];
        measure q[1];
        measure q[2];""", result_probability=(0.25, 0.25, 0, 0, 0.25, 0, 0, 0.25))  # 000, 001, 100, 111 all 0.25
    program_toffoli_with_flips = Program("""x q[0];
        x q[1];
        id q[2];
        h q[2];
        cx q[1], q[2];
        tdg q[2];
        cx q[0], q[2];
        t q[2];
        cx q[1], q[2];
        tdg q[2];
        cx q[0], q[2];
        t q[1];
        t q[2];
        h q[2];
        cx q[1], q[2];
        h q[1];
        h q[2];
        cx q[1], q[2];
        h q[1];
        h q[2];
        cx q[1], q[2];
        cx q[0], q[2];
        t q[0];
        tdg q[2];
        cx q[0], q[2];
        measure q[0];
        measure q[1];
        measure q[2];""", result_probability=(0, 0, 0, 0, 0, 0, 0, 1.0))  # 111: 1.0
    all_multi_gate_tests = [program_reverse_cnot, program_swap, program_swap_q0_q1, program_controlled_hadamard,
                            program_approximate_sqrtT, program_toffoli_state, program_toffoli_with_flips]
    # IBM Section IV, page 3 Grover's algorithm
    program_grover_n2_a00 = Program("""h q[1];
        h q[2];
        s q[1];
        s q[2];
        h q[2];
        cx q[1], q[2];
        h q[2];
        s q[1];
        s q[2];
        h q[1];
        h q[2];
        x q[1];
        x q[2];
        h q[2];
        cx q[1], q[2];
        h q[2];
        x q[1];
        x q[2];
        h q[1];
        h q[2];
        measure q[1];
        measure q[2];""", result_probability=(1.0, 0, 0, 0))  # 00: 1.0
    program_grover_n2_a01 = Program("""h q[1];
        h q[2];
        s q[2];
        h q[2];
        cx q[1], q[2];
        h q[2];
        s q[2];
        h q[1];
        h q[2];
        x q[1];
        x q[2];
        h q[2];
        cx q[1], q[2];
        h q[2];
        x q[1];
        x q[2];
        h q[1];
        h q[2];
        measure q[1];
        measure q[2];""", result_probability=(0.0, 1.0, 0.0, 0.0))  # 01: 1.0
    program_grover_n2_a10 = Program("""h q[1];
        h q[2];
        s q[1];
        h q[2];
        cx q[1], q[2];
        h q[2];
        s q[1];
        h q[1];
        h q[2];
        x q[1];
        x q[2];
        h q[2];
        cx q[1], q[2];
        h q[2];
        x q[1];
        x q[2];
        h q[1];
        h q[2];
        measure q[1];
        measure q[2];""", result_probability=(0, 0, 1.0, 0))  # 10: 1.0
    program_grover_n2_a11 = Program("""h q[1];
        h q[2];
        h q[2];
        cx q[1], q[2];
        h q[2];
        h q[1];
        h q[2];
        x q[1];
        x q[2];
        h q[2];
        cx q[1], q[2];
        h q[2];
        x q[1];
        x q[2];
        h q[1];
        h q[2];
        measure q[1];
        measure q[2];""", result_probability=(0, 0, 0, 1.0))  # 10: 1.0
    all_grover_tests = [program_grover_n2_a00, program_grover_n2_a01, program_grover_n2_a10, program_grover_n2_a11]
    # IBM Section IV, page 4 Deutsch-Jozsa Algorithm
    program_deutschjozsa_n3 = Program("""h q[0];
        h q[1];
        h q[2];
        h q[2];
        z q[0];
        cx q[1], q[2];
        h q[2];
        h q[0];
        h q[1];
        h q[2];
        measure q[0];
        measure q[1];
        measure q[2];""", result_probability=(0.25, 0.25, 0.25, 0.25))
    program_deutschjozsa_constant_n3 = Program("""h q[0];
        h q[1];
        h q[2];
        h q[0];
        h q[1];
        h q[2];
        measure q[0];
        measure q[1];
        measure q[2];""", result_probability=(1.0, 0, 0, 0))

    # IBM Section V, page 2 Quantum Repetition Code
    program_encoder_into_bitflip_code = Program("""h q[2];
        t q[2];
        h q[2];
        h q[1];
        h q[2];
        h q[3];
        cx q[1], q[2];
        cx q[3], q[2];
        h q[1];
        h q[2];
        h q[3];
        measure q[1];
        measure q[2];
        measure q[3];""", result_probability=(0.854, 0, 0, 0, 0, 0, 0, 0.146))
    program_encoder_and_decoder_tomography = Program("""h q[2];
        h q[1];
        h q[2];
        h q[3];
        cx q[1], q[2];
        cx q[3], q[2];
        h q[1];
        h q[2];
        h q[3];
        id q[1];
        id q[2];
        id q[3];
        id q[1];
        id q[2];
        id q[3];
        id q[1];
        id q[2];
        id q[3];
        h q[1];
        h q[2];
        h q[3];
        cx q[3], q[2];
        cx q[1], q[2];
        h q[1];
        h q[3];
        cx q[3], q[2];
        tdg q[2];
        cx q[1], q[2];
        t q[2];
        cx q[3], q[2];
        tdg q[2];
        cx q[1], q[2];
        t q[2];
        h q[2];
        bloch q[2];""", bloch_vals=(None, None, (1, 0, 0), None, None))  # Bloch q2: (1,0,0)
    program_encoder_into_bitflip_code_parity_checks = Program("""h q[2];
        t q[2];
        h q[2];
        h q[0];
        h q[1];
        h q[2];
        cx q[1], q[2];
        cx q[0], q[2];
        h q[0];
        h q[1];
        h q[3];
        cx q[3], q[2];
        h q[2];
        h q[3];
        cx q[3], q[2];
        cx q[0], q[2];
        cx q[1], q[2];
        h q[2];
        h q[4];
        cx q[4], q[2];
        h q[2];
        h q[4];
        cx q[4], q[2];
        cx q[1], q[2];
        cx q[3], q[2];
        measure q[2];
        measure q[4];
        measure q[0];
        measure q[1];
        measure q[3];""", result_probability=(
        0.852, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.146, 0, 0, 0, 0,
        0))  # 00000: 0.854; 11010: 0.146
    # IBM Section V, page 3 Stabilizer measurements
    program_plaquette_z0000 = Program("""id q[0];
        id q[1];
        id q[3];
        id q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(1.0, 0))
    program_plaquette_z0001 = Program("""id q[0];
        id q[1];
        id q[3];
        x q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(0, 1.0))
    program_plaquette_z0010 = Program("""id q[0];
        id q[1];
        x q[3];
        id q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(0, 1.0))
    program_plaquette_z0011 = Program("""id q[0];
        id q[1];
        x q[3];
        x q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(1.0, 0))
    program_plaquette_z0100 = Program("""id q[0];
        x q[1];
        id q[3];
        id q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(0, 1.0))
    program_plaquette_z0101 = Program("""id q[0];
        x q[1];
        id q[3];
        x q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(1.0, 0))
    program_plaquette_z0110 = Program("""id q[0];
        x q[1];
        x q[3];
        id q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(1.0, 0))
    program_plaquette_z0111 = Program("""id q[0];
        x q[1];
        x q[3];
        x q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(0, 1.0))
    program_plaquette_z1000 = Program("""x q[0];
        id q[1];
        id q[3];
        id q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(0, 1.0))
    program_plaquette_z1001 = Program("""x q[0];
        id q[1];
        id q[3];
        x q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(1.0, 0))
    program_plaquette_z1010 = Program("""x q[0];
        id q[1];
        x q[3];
        id q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(1.0, 0))
    program_plaquette_z1011 = Program("""x q[0];
        id q[1];
        x q[3];
        x q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(0, 1.0))
    program_plaquette_z1100 = Program("""x q[0];
        x q[1];
        id q[3];
        id q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(1.0, 0))
    program_plaquette_z1101 = Program("""x q[0];
        x q[1];
        id q[3];
        x q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(0, 1.0))
    program_plaquette_z1110 = Program("""x q[0];
        x q[1];
        x q[3];
        id q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(0, 1.0))
    program_plaquette_z1111 = Program("""x q[0];
        x q[1];
        x q[3];
        x q[4];
        cx q[4], q[2];
        cx q[0], q[2];
        cx q[3], q[2];
        cx q[1], q[2];
        measure q[2];""", result_probability=(1.0, 0))
    program_plaquette_zXplusminusplusminus = Program("""h q[0];
        h q[1];
        h q[3];
        h q[4];
        z q[1];
        z q[4];
        id q[0];
        id q[1];
        id q[3];
        id q[4];
        h q[0];
        h q[1];
        h q[3];
        h q[4];
        cx q[4], q[2];
        cx q[3], q[2];
        cx q[0], q[2];
        cx q[1], q[2];
        h q[0];
        h q[1];
        measure q[2];
        h q[3];
        h q[4];""", result_probability=(1.0, 0))
    # Convenience for testing
    all_normal_plaquette_programs = [program_plaquette_z0000, program_plaquette_z0001, program_plaquette_z0010,
                                     program_plaquette_z0011, program_plaquette_z0100, program_plaquette_z0101,
                                     program_plaquette_z0110, program_plaquette_z0111, program_plaquette_z1000,
                                     program_plaquette_z1001, program_plaquette_z1010, program_plaquette_z1011,
                                     program_plaquette_z1100, program_plaquette_z1101, program_plaquette_z1110,
                                     program_plaquette_z1111]
