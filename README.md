# REE-ED-Model

<h2>Introduction</h2>
<p>Mechanistic and probabilistic implementations of a model for the binary separation of REEs during electrodialysis. The code for 3 different implementations are included in this repository. The main implementation following the Takahashi et al. model is written in Python and includes a GUI, while an identical MATLAB implementation is also provided. To overcome the large discrepencies in experimental and model results, a machine learning model based off in-lab data was also created.</p>

<h2>Technologies</h2>
<p>The main Takahashi et al implementation requires a Jupyter notebook environment to be able to take advantage of the GUI. The code still runs on pure Python environments, with less flexibility in changing operating conditions after running the code. The machine learning model takes advantage of the TensorFlow library in Python to create a seamless model. The MATLAB implementation takes use of similar solving algorithims as the equivalent Python implementation.</p>

<h2>Examples of use</h2>
<p>For the Takahashi implementation, running the first 4 blocks of code will render a GUI with the following layout:</p>

![alt text](https://github.com/KellieChong/REE-ED-Model/blob/main/REE-ED-Takahashi-GUI.png)

<p>Once the GUI is run once, you can select the type of Rare Earth Elements in solution, as well as the chelating agent, and their respective concentrations. You can also select the pH, ion-exchange capacity of the membranes, number of time points to be iterated over, as well as the solution timespan using the sliders, or by directly enterring the value by double-clicking the existing value to access the text field. The solution will automatically update on the graphs each time an operating parameter is changed. The solver can also be selected by clicking one of the 3 numerical method solvers on the right. </p>
