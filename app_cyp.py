import gradio as gr

# Placeholder perfume database (can be replaced with real data later)
perfume_db = {
    "floral": ["Chanel No. 5", "Daisy by Marc Jacobs", "Flowerbomb by Viktor&Rolf"],
    "woody": ["Terre d’Hermès", "Tom Ford Oud Wood", "Santal 33"],
    "fresh": ["Light Blue by Dolce & Gabbana", "Cool Water", "CK One"],
    "spicy": ["YSL Opium", "Spicebomb", "Black Orchid"],
    "citrus": ["Clinique Happy", "Acqua di Parma", "Versace Man Eau Fraiche"]
}

# Logic: Based on user's preferences, recommend perfumes
def recommend_perfumes(preferences):
    recommendations = []
    for note in preferences:
        recommendations.extend(perfume_db.get(note, []))
    # Limit number of recommendations
    return list(set(recommendations))[:5] if recommendations else ["No matches found. Try selecting different notes."]

# Gradio UI with custom images and CSS in HTML/Markdown
with gr.Blocks() as demo:
    gr.Markdown("""
    <h1 style="text-align: center; color: #4a90e2;">🧴 Find Your Perfect Perfume</h1>
    <p style="text-align: center; color: #777;">Answer the following to get perfume recommendations:</p>
    """)
    
    fragrance_choices = gr.CheckboxGroup(
        ["floral", "woody", "fresh", "spicy", "citrus"],
        label="Which fragrance families do you enjoy?",
        info="You can select multiple.",
        elem_id="fragrance-selection"
    )


    gr.Markdown("""
    <style>
        #fragrance-selection {
            background-color: #f7f7f7;
            border-radius: 8px;
            padding: 10px;
            font-size: 16px;
        }
        
        #get-recommendations-btn {
            background-color: #4a90e2;
            color: white;
            font-size: 18px;
            padding: 10px 20px;
            border-radius: 5px;
            cursor: pointer;
            border: none;
            margin-top: 20px;
        }
        
        #get-recommendations-btn:hover {
            background-color: #357ABD;
        }

        #recommendations {
            background-color: #f7f7f7;
            border-radius: 8px;
            padding: 15px;
            margin-top: 20px;
            font-size: 16px;
            border: 1px solid #ddd;
        }

        h1 {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #4a90e2;
        }

        p {
            font-family: 'Arial', sans-serif;
            color: #555;
        }
    </style>
    """)

    output = gr.Textbox(label="Recommended Perfumes", elem_id="recommendations")
    get_button = gr.Button("Get Recommendations", elem_id="get-recommendations-btn")
    get_button.click(fn=recommend_perfumes, inputs=fragrance_choices, outputs=output)

demo.launch(share=False, inbrowser=True)

