import tkinter as tk
from tkinter import messagebox, colorchooser


# all what we need to handle preferences

class Preference:
    def __init__(self, name, value_type, default_value, string_conversion_method, description=""):
        self.name = name
        self.value_type = value_type
        self.default_value = default_value
        self.value = default_value
        self.string_conversion_method = string_conversion_method
        self.description = description
        self.preferencesPopup = None

    def set_value(self, new_value):
        try:
            # Only apply the string conversion method if the value is a string (this avoids issues with bools)
            if isinstance(new_value, str):
                self.value = self.string_conversion_method(new_value)
            else:
                self.value = new_value
        except ValueError:
            messagebox.showerror("Invalid Input", f"Invalid value for preference: {self.name}")
    
    def reset_to_default(self):
        self.value = self.default_value

    def get_value(self):
        return str(self.value)

        
    def dump(self):
        """Generates a detailed string dump of the preference."""
        dump = (
            f"Preference Dump:\n"
            f"Name: {self.name}\n"
            f"Name Type: {self.value_type.__name__}\n"
            f"Default Value: {self.default_value}\n"
            f"Current Value: {self.value}\n"
            f"Description: {self.description}\n"
        )
        print (dump)
        return dump   

class Preferences:
    def __init__(self, root, updateCallback:None):
        self.root = root
        self.preferences = self.initPreferences()
        self.updateCallback=updateCallback

    def closePreferencesPopup(self):
        self.preferencesPopup.destroy()

    def open_preferences_popup(self):
        def on_select(event):
            selected_index = listbox.curselection()
            prefKeysAsList=list(self.preferences)
            if selected_index:
                selectedPref:Preference = self.preferences.get(prefKeysAsList[selected_index[0]])
                print("name:",selectedPref.name," end value")
                # If the selected preference is a nucleotide color (A, T, C, G), directly open the color chooser
                if selectedPref.name in ['A', 'T', 'C', 'G']: 
                    self.pick_color(selectedPref, listbox)
                    self.updateCallback()   
                else:
                    self.edit_non_color_preference(selectedPref, listbox)
        self.preferencesPopup = tk.Toplevel(self.root)
        self.preferencesPopup.title("Preferences")
        # Make the preferences window modal
        self.preferencesPopup.grab_set()
        listbox = tk.Listbox(self.preferencesPopup, height=len(self.preferences), width=80, border=1)
        listbox.pack(padx=10, pady=10)
        self.update_preferences_display(listbox)
        listbox.bind("<Double-1>", on_select)
        close_button = tk.Button(self.preferencesPopup, text="Close", command=self.closePreferencesPopup)
        close_button.pack(pady=5)


    def update_preferences_display(self, listbox):
        """ Refresh the listbox display of preferences inside the popup. """
        listbox.delete(0, tk.END)
        for key, pref in self.preferences.items():
            if pref.value_type == bool:
                value = "Enabled" if pref.value else "Disabled"
                listbox.insert(tk.END, f"name: {pref.name}, value: {value}, description: {pref.description}, default: {pref.default_value}")                               
            else:
                listbox.insert(tk.END, f"name: {pref.name}, value: {pref.value}, description: {pref.description}, default: {pref.default_value}")

    def pick_color(self, preference, listbox):
        """ Open the color picker and update the color preference value. """
        color_code = colorchooser.askcolor(title=f"Choose Color for {preference.name}")[1]  # Get color in HEX format
        if color_code:
            self.update_preference_value(preference, color_code, listbox)

    def update_preference_value(self, preference, new_value, listbox):
        """ Update the preference value and refresh the UI. """
        preference.set_value(new_value)
        self.update_preferences_display(listbox)  # Update the listbox in the popup

    def edit_non_color_preference(self, preference:Preference, listbox):
        def apply_changes(event=None):
            if preference.value_type == bool:
                new_value = new_value_var.get()
                self.update_preference_value(preference, new_value, listbox)
            else:
                new_value = entry.get()
                self.update_preference_value(preference, new_value, listbox)

            # Update the listbox in the popup and main window after saving
            self.update_preferences_display(listbox)  # Update the popup's listbox
            self.updateCallback()   
            # Do not destroy popup immediately; keep it open and modal
            popup.grab_set()  # Reapply grab_set to ensure it remains modal
            popup.destroy()
            self.preferencesPopup.grab_set() # it seems to have amnesia and forgets it is modal

        def reset_to_default():
            preference.reset_to_default()
            self.update_preference_value(preference, preference.get_value(), listbox)
            # Update the listbox in the popup and main window after reset
            self.update_preferences_display(listbox)  # Update the popup's listbox
            popup.grab_set()  # Reapply grab_set to ensure it remains modal            
            popup.destroy()
            self.preferencesPopup.grab_set()  # it seems to have amnesia and forgets it is modal

        def commandCloseEditPopup():
            print(self)
            apply_changes()
            # popup.grab_set()  # Reapply grab_set to ensure it remains modal TODO adrian tkinter.TclError: bad window path name ".!toplevel2"           
            popup.destroy()             
            self.preferencesPopup.grab_set()  # it seems to have amnesia and forgets it is modal

        popup = tk.Toplevel(self.root)
        popup.title(f"Edit {preference.name}")

        # Make the edit preference window modal
        popup.grab_set()
        if preference.value_type == bool:
            new_value_var = tk.BooleanVar(value=preference.value)
            check_button = tk.Checkbutton(popup, text=f"{preference.name}: {preference.description}", variable=new_value_var)
            check_button.pack(padx=10, pady=5)

            save_button = tk.Button(popup, text="OK", command=apply_changes)
            save_button.pack(pady=5)
        else:
            entry = tk.Entry(popup)
            entry.insert(0, preference.get_value())
            entry.pack(padx=10, pady=5)

            # Bind the Enter key to apply changes
            entry.bind("<Return>", apply_changes)
            # No Save button for non-bool preferences, only Reset to Default
            reset_button = tk.Button(popup, text="OK", command=commandCloseEditPopup)
            reset_button.pack(pady=5)
            reset_button = tk.Button(popup, text="Reset to Default", command=reset_to_default)
            reset_button.pack(pady=5)            

    # def get_preference_value(self, preference_name):
    #     # Find the preference by name and return its current value
    #     p:Preference = self.preferences.get(preference_name)
    #     return p.get_value()

    def get_preference_value(self, preference_name):
        # Find the preference by name and return its current value
        p:Preference = self.preferences.get(preference_name)
        if p:
            v= p.get_value()
            if v is None:
                messagebox.showerror("Preference Not Found", f"Preference with name '{preference_name}' not found.")
            # elif p.value_type == bool:
            #     return True if v == "True" else False
            elif p.value_type != str:
                print("value before conversion",v)
                return p.string_conversion_method( v)
            return v
        messagebox.showerror("Preference Not Found", f"Preference with name '{preference_name}' not found.")
        return None
    
    def set_preference_value(self, preference_name, value):
        # Find the preference by name and return its current value
        p:Preference = self.preferences.get(preference_name)
        p.set_value(value)

    def initPreferences(self)->dict[str,Preference]:            
        def convertToBool(value:str)->bool:
            if value.lower() == "true":
                return True
            elif value.lower() == "false":
                return False
            else:
                raise ValueError("Invalid boolean value")
        # 
        p:dict[str,Preference] = {
            "defaultTestFileValue":     Preference("defaultTestFileValue", str,"/test.embl", str,"Path to the default file"),
            "canvasHorizontalMargin":   Preference("canvasHorizontalMargin", int,  0, int,"Horizontal margin for the canvas"),
            "fontName":                 Preference("fontName", str, "Arial",  str,"Font type for the canvas"),
            "fontSize":                 Preference("fontSize", int, 10, int,"Font type for the canvas"),
            "horizontalPixelsMargin":   Preference("horizontalPixelsMargin", int,  2,int, "Horizontal pixel margin for layout"),
            "verticalPixelsMargin":     Preference("verticalPixelsMargin", int, 2,  int,"Vertical pixel margin for layout"),
            "verticalSequenceSpacing":  Preference("verticalSequenceSpacing", int,  10, int,"Vertical sequence spacing"),
            "coloredBases":             Preference("coloredBases", bool, True, convertToBool, "Enable or disable colored bases in the sequence"),
            "shrink":             Preference("shrink", bool, True, convertToBool, "Shrinks the sequence and keeps features"),
            "A":                        Preference( "A", str, "cyan",  str,"Color for Adenine (A)"),
            "T":                        Preference( "T", str, "gold2",  str,"Color for Thymine (T)"),
            "G":                        Preference( "G", str, "lime green", str,"Color for Guanine (G)"),
            "C":                        Preference( "C", str, "red",  str,"Color for Cytosine (C)")
        }
        # print("name:",p.get("fontName").name," end value")
        return p
    
    def dump(self):
        for key, p in self.preferences.items():
            p.dump()


def main():
    root = tk.Tk()
    root.title("Preferences App")

    preferences=Preferences(root)
    open_button = tk.Button(root, text="Open Preferences", command=preferences.open_preferences_popup)
    print(preferences.get_preference_value("defaultTestFileValue"))
    open_button.pack(pady=20)
    root.mainloop()


if __name__ == "__main__":
    main()
