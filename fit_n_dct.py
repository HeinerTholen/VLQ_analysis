import os, itertools
import varial.generators as gen
from ROOT import TCanvas, TF1, gStyle, gROOT
gROOT.SetBatch()
gStyle.SetOptFit(0111)


#################################################################### general ###
def canvas_n_save_it(wrp):
    c = TCanvas(wrp.name, wrp.name, 600, 600)
    wrp.histo.Draw()
    tf1 = TF1("f", wrp.f1_str, 0., 10000.)
    tf1.Draw("same")
    c.SaveAs("plot_"+wrp.name+".png")
    return wrp


def make_html(wrps):
    names = list(w.name for w in wrps)
    f1_strs = list(
        'TF1* f_%s = new TF1("%s", "%s", 0., 1000.)<br>' % (
            w.name, w.name, w.f1_str
        ) for w in wrps
    )
    lines = ["<html>", "<h2>Names:</h2>"]
    lines += names
    lines += ["<h2>Function strings (same order):</h2>"]
    lines += f1_strs
    lines += ["<h2>Plots:</h2>"]
    lines += list("%s:<br><img src='plot_%s.png'>" % (s, s) for s in names)
    lines += ["</html>"]
    lines = list(l + "<br>\n" for l in lines)
    with open("index.html", "w") as f:
        f.writelines(lines)


############################################################### fit-specific ###
def htsub_do_da_fit(wrp):
    if "TTbar" in wrp.name:
        wrp.histo.Fit("expo")
        f1 = wrp.histo.GetFunction("expo")
        wrp.f1_str = "exp(%f*x + %f)" % (
            f1.GetParameter(1),
            f1.GetParameter(0)
        )
    else:
        wrp.histo.Fit("gaus")
        f1 = wrp.histo.GetFunction("gaus")
        wrp.f1_str = "%f * exp(-0.5 * ((x - %f )/ %f )**2)" % (
            f1.GetParameter(0),
            f1.GetParameter(1),
            f1.GetParameter(2)
        )
    return wrp


def mhigg_do_da_fit(wrp):
    peak = TF1("my_peak", "gaus", 0., 10000.)
    wrp.histo.Fit(peak)
    peakpol1 = TF1("my_peak_pol", "my_peak + pol1(3)", 0., 10000.)
    wrp.histo.Fit(peakpol1)
    wrp.f1_str = "%f * exp(-0.5 * ((x - %f )/ %f )**2) + %f + %f*x" % (
        peakpol1.GetParameter(0),
        peakpol1.GetParameter(1),
        peakpol1.GetParameter(2),
        peakpol1.GetParameter(3),
        peakpol1.GetParameter(4)
    )
    # peak = TF1("my_peak", "[0] / ((x^2 - [1]^2)^2 + ([1]*[2])^2)", 0., 10000.)
    # wrp.histo.Fit(peak)
    # gauspol1 = TF1("my_peak_pol", "my_peak + pol1(3)", 0., 10000.)
    # wrp.histo.Fit(gauspol1)
    # wrp.f1_str = "%f / ((x^2 - %f^2)^2 + (%f*%f)^2) + %f + %f*x" % (
    #     gauspol1.GetParameter(0),
    #     gauspol1.GetParameter(1),
    #     gauspol1.GetParameter(1),
    #     gauspol1.GetParameter(2),
    #     gauspol1.GetParameter(3),
    #     gauspol1.GetParameter(4)
    # )
    return wrp


def do_da_fit(wrp):
    if "mHigg" == wrp.name[:5]:
        return mhigg_do_da_fit(wrp)
    else:
        return htsub_do_da_fit(wrp)


def fit_get_funky(wrps):
    wrps = gen.gen_trim(wrps, True, False)
    wrps = (do_da_fit(w) for w in wrps)
    wrps = (canvas_n_save_it(w) for w in wrps)
    make_html(list(wrps))


############################################################### dct-specific ###
def make_dct_str(n, n_fade_out, histo):
    # see http://www.fftw.org/fftw3_doc/What-FFTW-Really-Computes.html
    values = list(histo.GetBinContent(i+1) for i in range(n))
    for i in range(1, n_fade_out):
        values[n - i] *= float(i) / n_fade_out
    s = "%f/2. + " % values[0]                  # x_0
    x_str = "(x-%.1f)/%.1f" % (
        histo.GetXaxis().GetXmin(),
        histo.GetXaxis().GetXmax() - histo.GetXaxis().GetXmin()
    )
    s += " + ".join(                            # x_1 -> x_n
        "%f*cos(pi*%d*%s)" % (values[i], i, x_str) for i in range(1, n)
    )
    s = "max(" + s + ", 0.00001)/%d." % histo.GetNbinsX()    # norm (fftw)
    return s


def do_da_dct(wrp):
    dct = gen.op.copy(wrp)
    wrp.histo.FFT(dct.histo, "RE R2R_2 M")
    wrp.f1_str = make_dct_str(11, 7, dct.histo)
    return wrp


def dct_get_funky(wrps):
    wrps = gen.gen_trim(wrps, True, False)
    wrps = (do_da_dct(w) for w in wrps)
    wrps = (canvas_n_save_it(w) for w in wrps)
    make_html(list(wrps))


####################################################################### main ###
def do_da_funk():
    wrps = gen.dir_content()
    wrps = itertools.ifilter(
        lambda w: (w.name[-2:] == "00"
                   or w.name[-3:] == "QCD"
                   or w.name[-5:] == "TTbar"),
        wrps
    )
    wrps = list(gen.load(wrps))
    if not os.path.exists("dct"):
        os.mkdir("dct")
    if not os.path.exists("fit"):
        os.mkdir("fit")
    os.chdir("fit")
    fit_get_funky(wrps)
    os.chdir("../dct")
    dct_get_funky(wrps)


if __name__ == "__main__":
    do_da_funk()
