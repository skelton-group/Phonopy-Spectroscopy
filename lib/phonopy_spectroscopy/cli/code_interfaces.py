# -*- coding: utf-8 -*-


class CodeInterface:
    def __init__(self):
        pass

    def load_gamma_ph(
        self,
        cell_file,
        freqs_evecs_file,
        lws_file=None,
        lws_t=300.0,
        irreps_file=None,
        conv_trans=None,
    ):
        raise NotImplementedError(
            "load_gamma_ph must be implemented in derived classes."
        )

    def load_born_eps_hf(self, born_file, eps_hf=None):
        raise NotImplementedError(
            "load_born_eps_hf must be implementd in derived classes."
        )

    def save_fd_raman_structures(self, fd_calc, output_fmt=None):
        raise NotImplementedError(
            "save_fd_raman_structures must be implemented in derived "
            "classes."
        )

    def read_fd_raman_dielectrics(self, fd_calc, input_files, input_fmt=None):
        raise NotImplementedError(
            "read_fd_raman_dielectrics must be implemented in derived "
            "classes."
        )

    def get_default_cell_file(self):
        return None

    def get_default_freqs_evecs_file(self):
        return None

    def get_default_irreps_file(self):
        return None


class PhonopyInterface(CodeInterface):
    pass
