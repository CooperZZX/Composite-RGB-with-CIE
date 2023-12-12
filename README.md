# Composite-RGB-with-CIE
Compotite true color images from hyperspectral data using CIE1931 standard
The hyperspectral images of CE-3, CE-4, and CE-5 are composited into true color images according the CIE 1931 standard by the following the steps:
1.Denosing: Denosing can help to reduce noise in the spectra. There are 65536 pixels of each hyperspectral image. Each pixel is a spectrum range from 450 nm to 945 nm @ 5 nm for CE-3 and CE-4, and 480 nm to 950 nm @ 5 nm for CE-5.
2.Interpolation: The CIE 1931 Standard (https://cie.co.at) starts at 380 nm and ends at 780 nm. Extrapolation of the hyperspectral data is needed in the short wave direction.
3.Dead Pixels Correction: There are statistic and dynamic dead pixels. The statistic dead pixels can be recognized manually. The dynamic dead pixels should be removed with proper algorithm dynamically.
True Color Image Composite: The RGB image is generated after obtaining the tristimulus values and color space coordinates according the CIE 1931 standard.
