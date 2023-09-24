python -m pytest --cov=spectral_wave_data --cov-report term-missing --html=report.html .
lynx report.html
