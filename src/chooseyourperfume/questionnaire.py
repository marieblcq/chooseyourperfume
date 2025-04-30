def ask_preferences():
    print("Welcome to the perfume selector!")
    family = input("Which fragrance family do you prefer? (Floral, Woody, Gourmand, Fruity, etc.): ")
    note = input("Do you prefer top notes, heart notes, or base notes?: ")
    return family.strip().lower(), note.strip().lower()