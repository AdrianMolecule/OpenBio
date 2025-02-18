## OpenBio: Open Source DNA & Amino Acid Sequence Manipulation Tool

OpenBio is an open-source application designed to manipulate DNA and amino acid sequences, filling a gap in available tools for sequence analysis. Unlike proprietary platforms like **Benchling**, **SnapGene**, or **Geneious**, which limit functionality and customization, OpenBio provides an extensible and open solution. It was created out of necessity to qualitatively verify **LAMP primers**, a task that lacked the right tool at the time.

OpenBio is built upon the powerful **BioPython** library, allowing users to perform various sequence manipulation tasks, including simulating **PCR** (Polymerase Chain Reaction) and **partial LAMP** reactions. Future updates will add support for additional features, such as **Gibson Assembly** simulations.

### Key Features:

- **DNA and Amino Acid Sequence Manipulation**: Easily work with DNA and amino acid sequences.
- **Simulator Functionality**: Simulate **PCR**, **partial LAMP**, and other molecular biology processes.
- **User-Friendly Visualization**: Features like text rotation for 3’ to 5’ strands and visual representation of hydrogen bonds make the app accessible for both **design** and **educational** purposes.
- **Educational and Research Use**: Ideal for use in educational settings as well as qualitative primer verification.

### Known Limitations:
As with any developing project, there are some limitations:

- **Supported Formats**: Currently, the tool supports the **GenBank (gb)** format. Unreleased code for **EMBL** format is available.
- **Feature Visualization**: Overlapping features in sequences don’t always display correctly, and visual representation of strands with two loops is still in development.
- **Primer Annealing**: At present, only one primer can be annealed at a time.

If you encounter any issues or would like to request improvements on these limitations, please reach out to the development team.

### Installation and Usage:

1. **Windows Executable**: For simplicity, a **precompiled Windows executable** (`OpenBio.exe`) is available for download. This eliminates the need for complex setup processes.

   Download the executable from the [**dist directory**](https://github.com/AdrianMolecule/OpenBio/blob/master/dist/OpenBio.exe).

2. **Running the Application**:
   - After downloading, simply click the `OpenBio.exe` file to launch the application on your Windows computer.
   - It may take a few seconds to start, but once running, the app operates quickly and efficiently.

3. **Source Code**:
   - The source code for **OpenBio** is available on [GitHub](https://github.com/AdrianMolecule/OpenBio). The project was developed using **Visual SourceSafe (VSS)**, making it relatively easy to install and contribute to.

### Development and Future Plans:

The application is developed based on the principle of adding features as needs arise. New capabilities and optimizations will be added based on user feedback and future requirements. Notably, **Gibson Assembly** simulation is next on the feature roadmap.

### Contact and Contributions:
If you encounter issues, have suggestions for improvement, or wish to contribute to the project, please feel free to contact the developer through the GitHub repository or via email.
