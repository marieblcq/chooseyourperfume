import streamlit as st
from questionnaire import scent_questionnaire

def main():
    st.title("🌸 Perfume Explorer App")

    menu = ["Home", "Scent Questionnaire"]
    choice = st.sidebar.selectbox("Navigation", menu)

    if choice == "Home":
        st.write("Welcome to the Perfume Explorer!")
    elif choice == "Scent Questionnaire":
        scent_questionnaire()

if __name__ == "__main__":
    main()