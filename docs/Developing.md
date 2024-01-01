# Design ideas

1. A pipeline has many *units*. For now, the RNA-seq pipeline was split into three units: QC, align, quantify, and the WGS pipeline was split into other three units: QC, align, and calling. No matter how the units change, two things remain constant:
- Users have to generate index for softwares.
- The input for the first unit must be provided by the user using a meta file, while others' are dependent on the previous unit output.

2. The required elements for executing each unit were split into three parts:
- Input and output file path for each each sample. This part is handled by `MetaPipe` class.
- Software related information e.g. software for this unit, software bin path. This part is handled by `ConfigPipe` class.
- Specific commands. This part is handled by class inherent from `PP` class.

I have designed this system in response to the flexibility of the software for each unit. By avoiding writing codes for each pipeline directly, we can significantly reduce the overall amount of code required. For example, the WGS workflow consists of only three units, but if each unit has two software options to choose from, it would result in a total of eight possible workflows.

Furthermore, in the event of adding new software in the future, it would require the creation of four additional pipelines. As the number of units and software increases, this will lead to an exponential growth in code.

However, under the current design, you only need:
- Write a new internal function to handle the input and output format in `MetaPipe` class
- Write a new internal function for generating index for this software in `IndexPipe` class
- Write a new internal funciton according to the software usage in `PP*` class.
- Add some keys in `ConfigPipe` attributes.

3. The design of `SampleIO` and `StepIO` bother me a lot. On the one hand, I want to standarize the interface to transfer input and output path between different units; On the other hand, the form of data in each unit is different. For example, the QC unit need fastq file, which has nothing to do with the calling unit. <u> Of course one choice is designing suitable data class for each unit, however, considering that a data class often only has one or two attributes, establishing a series of data classes in this way obviously complicates this simple task. </u> Through weighing various factor, I make `SampleIO` class more like an abstract base class (ABC): only has one fixed attribute: the sample name, and other attributes should be dependent on the exact unit. Since `StepIO` is based on `SampleIO` to some extent, it is also very flexible.