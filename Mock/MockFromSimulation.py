class MockFromGalaxia:

    def __init__(self, data_path):
        self._DataPath = data_path

    def generate_para_file(self, para_path, index_list, value_list, fn="output"):
        self._ParaPath = para_path
        self._paras_keys = [
            'outputFile',
            'outputDir',
            'photoSys',
            'magcolorNames',
            'appMagLimits[0]',
            'appMagLimits[1]',
            'absMagLimits[0]',
            'absMagLimits[1]',
            'colorLimits[0]',
            'colorLimits[1]',
            'geometryOption',
            'longitude',
            'latitude',
            'surveyArea',
            'fSample',
            'popID',
            'warpFlareOn',
            'seed',
            'r_max',
            'starType',
            'photoError'
        ]
        self._para_values = [
            'CSST_L000_B000',
            './',
            'SDSS',
            'g,g-r',
            10,
            11,
            -1000,
            1000,
            -1000,
            1000,
            1,
            0,
            0,
            2,
            1.0,
            -1,
            1,
            17,
            1000,
            0,
            0
        ]
        n_str = 36
        f = open(self._ParaPath + fn, 'w')
        for i_index, value in zip(index_list, value_list):
            self._para_values[i_index] = value
        n_para = len(self._para_values)
        for i in range(n_para):
            n_s = len(self._paras_keys[i])
            str_ws = (n_str-n_s)*' '
            f.write(self._paras_keys[i] +str_ws+f"{value_list[i]}\n")
        f.close()

    def run_galaxia(self, parafile):
        import subprocess
        import sys
        ret = subprocess.call(f"galaxia -r {parafile}", shell=True)
        if ret != 0:
            sys.exit('Exec cmd %s error, return value: %s' % (f"galaxia -r {parafile}", str(ret)))

    def ebf_to_fits(self, inputfile, outputfile, ind=None):
        import ebf
        print(inputfile, outputfile)
