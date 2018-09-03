# Author: corbett@caltech.edu

import numpy as np
import unittest
import re
import random
import itertools
from functools import reduce
from math import sqrt, pi, e, log
import time

#
# Gates
#

class Gate(object):
    i_ = np.complex(0, 1)
    ## One qubit gates
    # Hadamard gate
    H = 1. / sqrt(2) * np.matrix('1 1; 1 -1')
    # Pauli gates
    X = np.matrix('0 1; 1 0')
    Y = np.matrix([[0, -i_], [i_, 0]])
    Z = np.matrix([[1, 0], [0, -1]])

    # Defined as part of the Bell state experiment
    W = 1 / sqrt(2) * (X + Z)
    V = 1 / sqrt(2) * (-X + Z)

    # Other useful gates
    eye = np.eye(2, 2)

    S = np.matrix([[1, 0], [0, i_]])
    Sdagger = np.matrix([[1, 0], [0, -i_]])  # convenience Sdagger = S.conjugate().transpose()

    T = np.matrix([[1, 0], [0, e ** (i_ * pi / 4.)]])
    Tdagger = np.matrix([[1, 0], [0, e ** (-i_ * pi / 4.)]])  # convenience Tdagger= T.conjugate().transpose()

    # TODO: for CNOT gates define programatically instead of the more manual way below
    ## Two qubit gates
    # CNOT Gate (control is qubit 0, target qubit 1), this is the default CNOT gate
    CNOT2_01 = np.matrix('1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0')
    # control is qubit 1 target is qubit 0
    CNOT2_10 = np.kron(H, H) * CNOT2_01 * np.kron(H, H)  # =np.matrix('1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0')

    # operates on 2 out of 3 entangled qubits, control is first subscript, target second
    CNOT3_01 = np.kron(CNOT2_01, eye)
    CNOT3_10 = np.kron(CNOT2_10, eye)
    CNOT3_12 = np.kron(eye, CNOT2_01)
    CNOT3_21 = np.kron(eye, CNOT2_10)
    CNOT3_02 = np.matrix(
        '1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 1 0')
    CNOT3_20 = np.matrix(
        '1 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0; 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 1; 0 0 0 0 1 0 0 0; 0 1 0 0 0 0 0 0; 0 0 0 0 0 0 1 0; 0 0 0 1 0 0 0 0')

    # operates on 2 out of 4 entangled qubits, control is first subscript, target second
    CNOT4_01 = np.kron(CNOT3_01, eye)
    CNOT4_10 = np.kron(CNOT3_10, eye)
    CNOT4_12 = np.kron(CNOT3_12, eye)
    CNOT4_21 = np.kron(CNOT3_21, eye)
    CNOT4_13 = np.kron(eye, CNOT3_02)
    CNOT4_31 = np.kron(eye, CNOT3_20)
    CNOT4_02 = np.kron(CNOT3_02, eye)
    CNOT4_20 = np.kron(CNOT3_20, eye)
    CNOT4_23 = np.kron(eye, CNOT3_12)
    CNOT4_32 = np.kron(eye, CNOT3_21)
    CNOT4_03 = np.eye(16, 16)
    CNOT4_03[np.array([8, 9])] = CNOT4_03[np.array([9, 8])]
    CNOT4_03[np.array([10, 11])] = CNOT4_03[np.array([11, 10])]
    CNOT4_03[np.array([12, 13])] = CNOT4_03[np.array([13, 12])]
    CNOT4_03[np.array([14, 15])] = CNOT4_03[np.array([15, 14])]
    CNOT4_30 = np.eye(16, 16)
    CNOT4_30[np.array([1, 9])] = CNOT4_30[np.array([9, 1])]
    CNOT4_30[np.array([3, 11])] = CNOT4_30[np.array([11, 3])]
    CNOT4_30[np.array([5, 13])] = CNOT4_30[np.array([13, 5])]
    CNOT4_30[np.array([7, 15])] = CNOT4_30[np.array([15, 7])]

    # operates on 2 out of 5 entangled qubits, control is first subscript, target second
    CNOT5_01 = np.kron(CNOT4_01, eye)
    CNOT5_10 = np.kron(CNOT4_10, eye)
    CNOT5_02 = np.kron(CNOT4_02, eye)
    CNOT5_20 = np.kron(CNOT4_20, eye)
    CNOT5_03 = np.kron(CNOT4_03, eye)
    CNOT5_30 = np.kron(CNOT4_30, eye)
    CNOT5_12 = np.kron(CNOT4_12, eye)
    CNOT5_21 = np.kron(CNOT4_21, eye)
    CNOT5_13 = np.kron(CNOT4_13, eye)
    CNOT5_31 = np.kron(CNOT4_31, eye)
    CNOT5_14 = np.kron(eye, CNOT4_03)
    CNOT5_41 = np.kron(eye, CNOT4_30)
    CNOT5_23 = np.kron(CNOT4_23, eye)
    CNOT5_32 = np.kron(CNOT4_32, eye)
    CNOT5_24 = np.kron(eye, CNOT4_13)
    CNOT5_42 = np.kron(eye, CNOT4_31)
    CNOT5_34 = np.kron(eye, CNOT4_23)
    CNOT5_43 = np.kron(eye, CNOT4_32)
    CNOT5_04 = np.eye(32, 32)
    CNOT5_04[np.array([16, 17])] = CNOT5_04[np.array([17, 16])]
    CNOT5_04[np.array([18, 19])] = CNOT5_04[np.array([19, 18])]
    CNOT5_04[np.array([20, 21])] = CNOT5_04[np.array([21, 20])]
    CNOT5_04[np.array([22, 23])] = CNOT5_04[np.array([23, 22])]
    CNOT5_04[np.array([24, 25])] = CNOT5_04[np.array([25, 24])]
    CNOT5_04[np.array([26, 27])] = CNOT5_04[np.array([27, 26])]
    CNOT5_04[np.array([28, 29])] = CNOT5_04[np.array([29, 28])]
    CNOT5_04[np.array([30, 31])] = CNOT5_04[np.array([31, 30])]
    CNOT5_40 = np.eye(32, 32)
    CNOT5_40[np.array([1, 17])] = CNOT5_40[np.array([17, 1])]
    CNOT5_40[np.array([3, 19])] = CNOT5_40[np.array([19, 3])]
    CNOT5_40[np.array([5, 21])] = CNOT5_40[np.array([21, 5])]
    CNOT5_40[np.array([7, 23])] = CNOT5_40[np.array([23, 7])]
    CNOT5_40[np.array([9, 25])] = CNOT5_40[np.array([25, 9])]
    CNOT5_40[np.array([11, 27])] = CNOT5_40[np.array([27, 11])]
    CNOT5_40[np.array([13, 29])] = CNOT5_40[np.array([29, 13])]
    CNOT5_40[np.array([15, 31])] = CNOT5_40[np.array([31, 15])]


#
# States
#

class State(object):
    i_ = np.complex(0, 1)
    ## One qubit states (basis)
    # standard basis (z)
    zero_state = np.matrix('1; 0')
    one_state = np.matrix('0; 1')
    # diagonal basis (x)
    plus_state = 1 / sqrt(2) * np.matrix('1; 1')
    minus_state = 1 / sqrt(2) * np.matrix('1; -1')
    # circular basis (y)
    plusi_state = 1 / sqrt(2) * np.matrix([[1], [i_]])  # also known as clockwise arrow state
    minusi_state = 1 / sqrt(2) * np.matrix([[1], [-i_]])  # also known as counterclockwise arrow state

    # 2-qubit states
    bell_state = 1 / sqrt(2) * np.matrix('1; 0; 0; 1')

    @staticmethod
    def change_to_x_basis(state):
        return Gate.H * state

    @staticmethod
    def change_to_y_basis(state):
        return Gate.H * Gate.Sdagger * state

    @staticmethod
    def change_to_w_basis(state):
        # W=1/sqrt(2)*(X+Z)
        return Gate.H * Gate.T * Gate.H * Gate.S * state

    @staticmethod
    def change_to_v_basis(state):
        # V=1/sqrt(2)*(-X+Z)
        return Gate.H * Gate.Tdagger * Gate.H * Gate.S * state

    @staticmethod
    def is_fully_separable(qubit_state):
        try:
            separated_state = State.separate_state(qubit_state)
            for state in separated_state:
                State.string_from_state(state)
            return True
        except StateNotSeparableException as e:
            return False

    @staticmethod
    def get_first_qubit(qubit_state):
        return State.separate_state(qubit_state)[0]

    @staticmethod
    def get_second_qubit(qubit_state):
        return State.separate_state(qubit_state)[1]

    @staticmethod
    def get_third_qubit(qubit_state):
        return State.separate_state(qubit_state)[2]

    @staticmethod
    def get_fourth_qubit(qubit_state):
        return State.separate_state(qubit_state)[3]

    @staticmethod
    def get_fifth_qubit(qubit_state):
        return State.separate_state(qubit_state)[4]

    @staticmethod
    def all_state_strings(n_qubits):
        return [''.join(map(str, state_desc)) for state_desc in itertools.product([0, 1], repeat=n_qubits)]

    @staticmethod
    def state_from_string(qubit_state_string):
        if not all(x in '01' for x in qubit_state_string):
            raise Exception("Description must be a string in binary")
        state = None
        for qubit in qubit_state_string:
            if qubit == '0':
                new_contrib = State.zero_state
            elif qubit == '1':
                new_contrib = State.one_state
            if state is None:
                state = new_contrib
            else:
                state = np.kron(state, new_contrib)
        return state

    @staticmethod
    def string_from_state(qubit_state):
        separated = State.separate_state(qubit_state)
        desc = ''
        for state in separated:
            if np.allclose(state, State.zero_state):
                desc += '0'
            elif np.allclose(state, State.one_state):
                desc += '1'
            else:
                raise StateNotSeparableException("State is not separable")
        return desc

    @staticmethod
    def separate_state(qubit_state):
        """This only works if the state is fully separable at present

        Throws exception if not a separable state"""
        n_entangled = QuantumRegister.num_qubits(qubit_state)
        if list(qubit_state.flat).count(1) == 1:
            separated_state = []
            idx_state = list(qubit_state.flat).index(1)
            add_factor = 0
            size = qubit_state.shape[0]
            while (len(separated_state) < n_entangled):
                size = size / 2
                if idx_state < (add_factor + size):
                    separated_state += [State.zero_state]
                    add_factor += 0
                else:
                    separated_state += [State.one_state]
                    add_factor += size
            return separated_state
        else:
            # Try a few naive separations before giving up
            cardinal_states = [State.zero_state, State.one_state, State.plus_state, State.minus_state,
                               State.plusi_state, State.minusi_state]
            for separated_state in itertools.product(cardinal_states, repeat=n_entangled):
                candidate_state = reduce(lambda x, y: np.kron(x, y), separated_state)
                if np.allclose(candidate_state, qubit_state):
                    return separated_state
            # TODO: more general separation methods
            raise StateNotSeparableException(
                "TODO: Entangled qubits not represented yet in quantum computer implementation. Can currently do manual calculations; see TestBellState for Examples")

    @staticmethod
    def measure(state):
        """finally some probabilities, whee. To properly use, set the qubit you measure to the result of this function
            to collapse it. state=measure(state). Currently supports only up to three entangled qubits """
        state_z = state
        n_qubits = QuantumRegister.num_qubits(state)
        probs = Probability.get_probabilities(state_z)
        rand = random.random()
        for idx, state_desc in enumerate(State.all_state_strings(n_qubits)):
            if rand < sum(probs[0:(idx + 1)]):
                return State.state_from_string(state_desc)

    @staticmethod
    def get_bloch(state):
        return np.array(
            (Probability.expectation_x(state), Probability.expectation_y(state), Probability.expectation_z(state)))

    @staticmethod
    def pretty_print_gate_action(gate, n_qubits):
        for s in list(itertools.product([0, 1], repeat=n_qubits)):
            sname = ('%d' * n_qubits) % s
            state = State.state_from_string(sname)
            print(sname, '->', State.string_from_state(gate * state))


class StateNotSeparableException(Exception):
    def __init__(self, args=None):
        self.args = args


class Probability(object):
    @staticmethod
    def get_probability(coeff):
        return (coeff * coeff.conjugate()).real

    @staticmethod
    def get_probabilities(state):
        return [Probability.get_probability(x) for x in state.flat]

    @staticmethod
    def get_correlated_expectation(state):
        probs = Probability.get_probabilities(state)
        return probs[0] + probs[3] - probs[1] - probs[2]

    @staticmethod
    def pretty_print_probabilities(state):
        probs = Probability.get_probabilities(state)
        am_desc = '|psi>='
        pr_desc = ''
        for am, pr, state_desc in zip(state.flat, probs, State.all_state_strings(QuantumRegister.num_qubits(state))):
            if am != 0:
                if am != 1:
                    am_desc += '%r|%s>+' % (am, state_desc)
                else:
                    am_desc += '|%s>+' % (state_desc)
            if pr != 0:
                pr_desc += 'Pr(|%s>)=%f; ' % (state_desc, pr)
        print(am_desc[0:-1])
        print(pr_desc)
        if state.shape == (4, 1):
            print("<state>=%f" % float(probs[0] + probs[3] - probs[1] - probs[2]))

    @staticmethod
    def expectation_x(state):
        state_x = State.change_to_x_basis(state)
        prob_zero_state = (state_x.item(0) * state_x.item(0).conjugate()).real
        prob_one_state = (state_x.item(1) * state_x.item(1).conjugate()).real
        return prob_zero_state - prob_one_state

    @staticmethod
    def expectation_y(state):
        state_y = State.change_to_y_basis(state)
        prob_zero_state = (state_y.item(0) * state_y.item(0).conjugate()).real
        prob_one_state = (state_y.item(1) * state_y.item(1).conjugate()).real
        return prob_zero_state - prob_one_state

    @staticmethod
    def expectation_z(state):
        state_z = state
        prob_zero_state = (state_z.item(0) * state_z.item(0).conjugate()).real
        prob_one_state = (state_z.item(1) * state_z.item(1).conjugate()).real
        return prob_zero_state - prob_one_state


class QuantumRegister(object):
    def __init__(self, name, state=State.zero_state, entangled=None):
        self._entangled = [self]
        self._state = state
        self.name = name
        self.idx = None
        self._noop = []  # after a measurement set this so that we can allow no further operations. Set to Bloch coords if bloch operation performed

    @staticmethod
    def num_qubits(state):
        num_qubits = log(state.shape[0], 2)
        if state.shape[1] != 1 or num_qubits not in [1, 2, 3, 4, 5]:
            raise Exception("unrecognized state shape")
        else:
            return int(num_qubits)

    def get_entangled(self):
        return self._entangled

    def set_entangled(self, entangled):
        self._entangled = entangled
        for qb in self._entangled:
            qb._state = self._state
            qb._entangled = self._entangled

    def get_state(self):
        return self._state

    def set_state(self, state):
        self._state = state
        for qb in self._entangled:
            qb._state = state
            qb._entangled = self._entangled
            qb._noop = self._noop

    def get_noop(self):
        return self._noop

    def set_noop(self, noop):
        self._noop = noop
        for qb in self._entangled:
            qb._noop = noop

    def is_entangled(self):
        return len(self._entangled) > 1

    def is_entangled_with(self, qubit):
        return qubit in self._entangled

    def get_indices(self, target_qubit):
        if not self.is_entangled_with(target_qubit):
            search = self._entangled + target_qubit.get_entangled()
        else:
            search = self._entangled
        return search.index(self), search.index(target_qubit)

    def get_num_qubits(self):
        return QuantumRegister.num_qubits(self._state)

    def __eq__(self, other):
        if not isinstance(other, type(self)): return NotImplemented
        return self.name == other.name and np.array(self._noop).shape == np.array(other._noop).shape and np.allclose(
            self._noop, other._noop) and np.array(self.get_state()).shape == np.array(
            other.get_state()).shape and np.allclose(self.get_state(),
                                                     other.get_state()) and QuantumRegisterCollection.orderings_equal(
            self._entangled, other._entangled)


class QuantumRegisterSet(object):
    """Created this so I could have some set like features for use, even though QuantumRegisters are mutable"""
    registers = []

    def __init__(self, registers):
        for r in registers:
            if r not in self.registers:
                self.registers += [r]

    def intersection(self, quantumregisterset):
        intersection = []

        if self.size() >= quantumregisterset:
            qrs1 = self
            qrs2 = quantumregisterset
        else:
            qrs1 = quantumregisterset
            qrs2 = self
        # now qrs2 is the smaller set
        intersection = [qr for qr in qrs1 if qr in qrs2]
        return QuantumRegisterSet(intersection)

    def size(self):
        return len(self.registers)


class QuantumRegisterCollection(object):
    def __init__(self, qubits):
        self._qubits = qubits
        for idx, qb in enumerate(self._qubits):
            qb.idx = idx
        self.num_qubits = len(qubits)

    def get_quantum_register_containing(self, name):
        for qb in self._qubits:
            if qb.name == name:
                return qb
            else:
                for entqb in qb.get_entangled():
                    if entqb.name == name:
                        return entqb
        raise Exception("qubit %s not found in %s" % (name, repr(self._qubits)))

    def get_quantum_registers(self):
        return self._qubits

    def entangle_quantum_registers(self, first_qubit, second_qubit):
        new_entangle = first_qubit.get_entangled() + second_qubit.get_entangled()
        if len(first_qubit.get_entangled()) >= len(second_qubit.get_entangled()):
            self._remove_quantum_register_named(second_qubit.name)
            first_qubit.set_entangled(new_entangle)
        else:
            self._remove_quantum_register_named(first_qubit.name)
            second_qubit.set_entangled(new_entangle)

    def _remove_quantum_register_named(self, name):
        self._qubits = [qb for qb in self._qubits if qb.name != name]

    def is_in_canonical_ordering(self):
        return self.get_qubit_order() == list(range(self.num_qubits))

    @staticmethod
    def is_in_increasing_order(qb_list):
        for a, b in zip(qb_list, qb_list[1:]):
            if not a.idx < b.idx:
                return False
        return True

    def get_entangled_qubit_order(self):
        ordering = []
        for qb in self._qubits:
            ent_order = []
            for ent in qb.get_entangled():
                ent_order += [ent]
            ordering += [ent_order]
        return ordering

    def get_qubit_order(self):
        ordering = []
        for qb in self._qubits:
            for ent in qb.get_entangled():
                ordering += [ent.idx]
        return ordering

    def add_quantum_register(self, qubit):
        qubit.idx = self.num_qubits
        self._qubits += [qubit]
        self.num_qubits += 1

    @staticmethod
    def orderings_equal(order_one, order_two):
        return [qb.idx for qb in order_one] == [qb.idx for qb in order_two]


class QuantumComputer(object):
    """This class is meant to simulate the 5-qubit IBM quantum computer,
        and be able to interpret the auto generated programs on the site.

        For entangled states, qubits are always reported in alphanumerical order
        """

    def __init__(self):
        self.qubits = QuantumRegisterCollection(
            [QuantumRegister("q0"), QuantumRegister("q1"), QuantumRegister("q2"), QuantumRegister("q3"),
             QuantumRegister("q4")])

    def reset(self):
        self.qubits = QuantumRegisterCollection(
            [QuantumRegister("q0"), QuantumRegister("q1"), QuantumRegister("q2"), QuantumRegister("q3"),
             QuantumRegister("q4")])

    def get_ordering(self):
        return self.qubits.get_qubit_order()

    def is_in_canonical_ordering(self):
        return self.qubits.is_in_canonical_ordering()

    def get_requested_state_order(self, name):
        get_states_for = [self.qubits.get_quantum_register_containing(x.strip()) for x in name.split(',')]
        if not QuantumRegisterCollection.is_in_increasing_order(get_states_for):
            raise Exception("at this time, requested qubits must be in increasing order")
        entangled_qubit_order = self.qubits.get_entangled_qubit_order()
        # # We know the idxs run range(5)
        # # We know if the idxs are contiguous, increasing we are good
        for get_state_for_qb in get_states_for:
            for eqb in entangled_qubit_order:
                eqo = [q.idx for q in eqb]
                # We know if the idxs are missing a number AND we want to find an idx that lies in there, we must entangle those states
                if not get_state_for_qb.idx in eqo and get_state_for_qb.idx in range(min(eqo), max(eqo) + 1):
                    print("We'll have to entangle the two")
                    # We'll have to entangle the two
                    qb1 = self.qubits.get_quantum_register_containing(eqo[0].name)
                    get_state_for_qb.set_state(np.kron(qb.get_state(), qb1.get_state()))
                    self.qubits.entangle_quantum_registers(get_state_for_qb, qb1)
                    return self.qubit_states_equal(name, state)

        # OK, if we reach here, we have all the entanglement we need, and we just need to sort the individual entangled states to match the output order
        for qubit in self.qubits.get_quantum_registers():
            if not QuantumRegisterCollection.is_in_increasing_order(qubit.get_entangled()):  # all one apart
                # We're not in order
                # We need to assert that the full return can be comprised of concatenating states from beginning to end without extras
                if not QuantumRegisterSet(qubit.get_entangled()).size() <= QuantumRegisterSet(
                        get_states_for).size() and QuantumRegisterSet(qubit.get_entangled()).intersection(
                        QuantumRegisterSet(get_states_for)).size():
                    raise Exception(
                        "With this entanglement setup we can't fully separate out just the qubits of iterest. Try measuring more bits")
                # We only care if we actually want to return something from this state Put eqo in order then
                # We want a sorting algorithm that easily maps to matrix operations, since we only have 5 elements max
                # we'll use bubble sort
                swapped = True
                n = len(qubit.get_entangled())
                while (swapped):
                    swapped = False
                    current_entangled = qubit.get_entangled()

                    for idx in range(len(current_entangled) - 1):
                        first_qubit = current_entangled[idx]
                        second_qubit = current_entangled[idx + 1]
                        if first_qubit.idx > second_qubit.idx:
                            current_entangled[idx] = second_qubit
                            current_entangled[idx + 1] = first_qubit
                            permute = np.eye(2 ** n, 2 ** n)
                            all_combos = list(itertools.product([0, 1], repeat=n))
                            already_swapped = []
                            for icombo, combo in enumerate(all_combos[:len(all_combos)]):
                                new_combo = list(combo)
                                new_combo[idx] = combo[idx + 1]
                                new_combo[idx + 1] = combo[idx]
                                new_combo = tuple(new_combo)
                                if combo != new_combo:
                                    inew_combo = all_combos.index(new_combo)
                                    swapset = set([icombo, inew_combo])
                                    if not swapset in already_swapped:
                                        already_swapped += [swapset]
                                        permute[np.array([icombo, inew_combo])] = permute[
                                            np.array([inew_combo, icombo])]
                            first_qubit.set_entangled(current_entangled)
                            first_qubit.set_state(permute * first_qubit.get_state())
                            swapped = True
        # OK, if we reach here, everything is in order, and entangled states are either all of interest or none are of interest we just need to return it!
        answer_state = None
        for qb in self.qubits.get_quantum_registers():
            if QuantumRegisterSet(qb.get_entangled()).size() <= QuantumRegisterSet(get_states_for).size():
                if answer_state is None:
                    answer_state = qb.get_state()
                else:
                    answer_state = np.kron(answer_state, qb.get_state())
        return answer_state

    def probabilities_equal(self, name, prob):
        get_states_for = [self.qubits.get_quantum_register_containing(x.strip()) for x in name.split(',')]
        if not QuantumRegisterCollection.is_in_increasing_order(get_states_for):
            raise Exception("at this time, requested qubits must be in increasing order")
        entangled_qubit_order = self.qubits.get_entangled_qubit_order()
        if (len(get_states_for) == 1 and self.is_in_canonical_ordering()) or (
                [x.name for x in get_states_for] in [[x.name for x in l] for l in entangled_qubit_order]):
            return np.allclose(Probability.get_probabilities(get_states_for[0].get_state()), prob)
        else:
            answer_state = self.get_requested_state_order(name)
            return np.allclose(Probability.get_probabilities(answer_state), prob, atol=1e-2)

    def qubit_states_equal(self, name, state):
        get_states_for = [self.qubits.get_quantum_register_containing(x.strip()) for x in name.split(',')]
        if not QuantumRegisterCollection.is_in_increasing_order(get_states_for):
            raise Exception("at this time, requested qubits must be in increasing order")
        entangled_qubit_order = self.qubits.get_entangled_qubit_order()
        if (len(get_states_for) == 1 and self.is_in_canonical_ordering()) or (get_states_for in entangled_qubit_order):
            return np.allclose(get_states_for[0].get_state(), state)
        else:
            answer_state = self.get_requested_state_order(name)
            return np.allclose(answer_state, state)

    def bloch_coords_equal(self, name, coords):
        on_qubit = self.qubits.get_quantum_register_containing(name)
        if self.is_in_canonical_ordering() and not on_qubit.is_entangled():
            return np.allclose(on_qubit.get_noop(), coords, atol=1e-3)
        else:
            try:
                separated_qubit = State.separate_state(on_qubit.get_state())
                on_qubit_idx = (on_qubit.get_entangled()).index(on_qubit)
                return np.allclose(State.get_bloch(separated_qubit[on_qubit_idx]), coords, atol=1e-3)
            except StateNotSeparableException as e:
                raise Exception("Entangled measurements that cannot be separatednot yet implemented for bloch sphere")

    def apply_gate(self, gate, on_qubit_name):
        on_qubit = self.qubits.get_quantum_register_containing(on_qubit_name)
        if len(on_qubit.get_noop()) > 0:
            print(
                "NOTE this qubit has been measured previously, there should be no more gates allowed but we are reverting that measurement for consistency with IBM's language")
            on_qubit.set_state(on_qubit.get_noop())
            on_qubit.set_noop([])
        if not on_qubit.is_entangled():
            if on_qubit.get_num_qubits() != 1:
                raise Exception("This qubit is not marked as entangled but it has an entangled state")
            on_qubit.set_state(gate * on_qubit.get_state())
        else:
            if not on_qubit.get_num_qubits() > 1:
                raise Exception("This qubit is marked as entangled but it does not have an entangled state")
            n_entangled = len(on_qubit.get_entangled())
            apply_gate_to_qubit_idx = [qb.name for qb in on_qubit.get_entangled()].index(on_qubit_name)
            if apply_gate_to_qubit_idx == 0:
                entangled_gate = gate
            else:
                entangled_gate = Gate.eye
            for i in range(1, n_entangled):
                if apply_gate_to_qubit_idx == i:
                    entangled_gate = np.kron(entangled_gate, gate)
                else:
                    entangled_gate = np.kron(entangled_gate, Gate.eye)
            on_qubit.set_state(entangled_gate * on_qubit.get_state())

    def apply_two_qubit_gate_CNOT(self, first_qubit_name, second_qubit_name):
        """ Should work for all combination of qubits"""
        first_qubit = self.qubits.get_quantum_register_containing(first_qubit_name)
        second_qubit = self.qubits.get_quantum_register_containing(second_qubit_name)
        if len(first_qubit.get_noop()) > 0 or len(second_qubit.get_noop()) > 0:
            raise Exception("Control or target qubit has been measured previously, no more gates allowed")
        if not first_qubit.is_entangled() and not second_qubit.is_entangled():
            combined_state = np.kron(first_qubit.get_state(), second_qubit.get_state())
            if first_qubit.get_num_qubits() != 1 or second_qubit.get_num_qubits() != 1:
                raise Exception("Both qubits are marked as not entangled but one or the other has an entangled state")
            new_state = Gate.CNOT2_01 * combined_state
            if State.is_fully_separable(new_state):
                second_qubit.set_state(State.get_second_qubit(new_state))
            else:
                self.qubits.entangle_quantum_registers(first_qubit, second_qubit)
                first_qubit.set_state(new_state)
        else:
            if not first_qubit.is_entangled_with(second_qubit):
                # Entangle the state
                combined_state = np.kron(first_qubit.get_state(), second_qubit.get_state())
                self.qubits.entangle_quantum_registers(first_qubit, second_qubit)
            else:
                # We are ready to do the operation
                combined_state = first_qubit.get_state()
            # Time for more meta programming!
            # Select gate based on indices
            control_qubit_idx, target_qubit_idx = first_qubit.get_indices(second_qubit)
            gate_size = QuantumRegister.num_qubits(combined_state)
            try:
                namespace = locals()
                exec('gate=Gate.CNOT%d_%d%d' % (gate_size, control_qubit_idx, target_qubit_idx), globals(), namespace)
                gate = namespace['gate']
            except:
                print('gate=Gate.CNOT%d_%d%d' % (gate_size, control_qubit_idx, target_qubit_idx))
                raise Exception("Unrecognized combination of number of qubits")
            first_qubit.set_state(gate * combined_state)

    def bloch(self, qubit_name):
        on_qubit = self.qubits.get_quantum_register_containing(qubit_name)
        if len(on_qubit.get_noop()) == 0:
            if not on_qubit.is_entangled():
                on_qubit.set_noop(State.get_bloch(on_qubit.get_state()))
            else:
                on_qubit.set_noop([1])

    def measure(self, qubit_name):
        on_qubit = self.qubits.get_quantum_register_containing(qubit_name)
        if len(on_qubit.get_noop()) == 0:
            on_qubit.set_noop(on_qubit.get_state())  # state before measurement for testing
            on_qubit.set_state(State.measure(on_qubit.get_state()))

    def execute(self, program):
        """Time for some very lazy meta programming!
        """
        # Transforming IBM's language to my variables
        lines = program.split(';')
        translation = [
            ['q[0]', '"q0"'],
            ['q[1]', '"q1"'],
            ['q[2]', '"q2"'],
            ['q[3]', '"q3"'],
            ['q[4]', '"q4"'],
            ['bloch ', r'self.bloch('],
            ['measure ', r'self.measure('],
            ['id ', 'self.apply_gate(Gate.eye,'],
            ['sdg ', 'self.apply_gate(Gate.Sdagger,'],
            ['tdg ', 'self.apply_gate(Gate.Tdagger,'],
            ['h ', 'self.apply_gate(Gate.H,'],
            ['t ', 'self.apply_gate(Gate.T,'],
            ['s ', 'self.apply_gate(Gate.S,'],
            ['x ', 'self.apply_gate(Gate.X,'],
            ['y ', 'self.apply_gate(Gate.Y,'],
            ['z ', 'self.apply_gate(Gate.Z,'],
        ]
        cnot_re = re.compile('^cx (q\[[0-4]\]), (q\[[0-4]\])$')
        for l in lines:
            l = l.strip()
            if not l: continue
            # CNOT operates on two qubits so gets special processing
            cnot = cnot_re.match(l)
            if cnot:
                control_qubit = cnot.group(1)
                target_qubit = cnot.group(2)
                l = 'self.apply_two_qubit_gate_CNOT(%s,%s' % (control_qubit, target_qubit)
            for k, v in translation:
                l = l.replace(k, v)
            l = l + ')'
            # Now running the code
            exec(l, globals(), locals())


class Program(object):
    def __init__(self, code, result_probability=[], bloch_vals=()):
        self.code = code
        self.result_probability = result_probability
        self.bloch_vals = bloch_vals
