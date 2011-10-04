

def add_help_text(form, help_text):
    for key,field in form.fields.items():
        field.help_text = help_text.get(key,'')
