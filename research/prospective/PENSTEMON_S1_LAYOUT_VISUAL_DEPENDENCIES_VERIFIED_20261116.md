# Original supplement source-render retry after dependency audit

The original PDF source loader imports `pypdf`, while the visual renderer requires PyMuPDF. The first one-shot job may fail before source inspection because only PyMuPDF had been listed. This retry installs both packages, retains the same original PDF SHA256, and makes no biological-source or eligibility rule change. The only output is temporary, original-source image evidence for the 52-versus-51 voucher discrepancy.
