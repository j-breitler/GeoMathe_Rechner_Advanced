# 📐 GeoMathe Rechner

Ein Python‑basierter Rechner für Aufgaben aus **GeoMathe 1** mit mehreren Rechenmodulen und zwei verfügbaren Versionen.

---

## 👤 Autor:innen

* **Hauptautor:innen:** Ai & Julian Breitler
* **Co‑Autor:innen & Unterstützung:** Fenja Runfors, Klara Gössler, Christopher Stering

---

## ⚠️ Haftungsausschluss

**Benützung auf eigene Gefahr!**
Für fehlerhafte Ergebnisse oder daraus entstehende Konsequenzen wird keine Haftung übernommen.

---

## 🛠 Voraussetzungen

* Installiertes **Python** (empfohlen: aktuelle Version)

---

## 📦 Versionen

### 🔹 Version 1 – *GeoMathe_Rechner*

* Verwendet ausschließlich die **Python‑Standardbibliothek**
* Keine zusätzlichen Packages notwendig
* Einfachste und sicherste Variante

### 🔹 Version 2 – *GeoMathe_Rechner_Advanced*

* Modernes und übersichtliches **User Interface**
  (größere Schrift, sauberes Layout)
* Kann eingebene Punkte zeichnen und in einem Koordinatensystem darstellen (sehr Hilfreich für Skizze!)
* Benötigt zusätzlich das externe Package **matplotlib**
* Installation des Packages erfolgt je nach System manuell oder automatisch

> 💡 Tipp: Wie man Python‑Packages installiert, kann einfach gegoogelt oder mithilfe einer KI nachgeschlagen werden.

---

## ▶️ Starten des Rechners

Beide Versionen können auf mehreren Wegen gestartet werden:

* 📂 **Doppelklick** auf die `.py`‑Datei
* 💻 Start über **Terminal / Command Line**
* 🧠 Ausführen in einer IDE (z. B. *VS Code, PyCharm, Spyder*)
* 🪟 **`.exe`‑Datei** per Doppelklick (nur Advanced‑Version, **Windows‑only**, installiert `matplotlib` automatisch)

### ✅ Empfehlung

Die einfachste und zuverlässigste Methode ist der **Doppelklick auf die `.py`‑Datei** der einfachen Version
(*funktioniert auf allen Systemen mit installiertem Python*).

---

## 🧭 Programmaufbau

Der Rechner ist in **fünf Tabs** gegliedert:

---

### 🟦 1. Tab – Point Management

* Eingabe von Punkten
* Speicherung unter einer **Point ID**
* Gespeicherte Punkte können in **Tab 2 (2.HA)** und **Tab 3 (1.HA)** wiederverwendet werden

---

### 🟩 2. Tab – 2. Hauptaufgabe

* Auswahl gespeicherter Punkte über ein **Dropdown‑Menü**
* Automatische **Modulo‑Operation**, falls:

  * Richtung > 400 gon
  * Richtung negativ ist
* Originalwert bleibt sichtbar, Hinweis zur Modulo‑Operation wird angezeigt

**Code‑Referenz:**

* *GeoMathe_Rechner:* Zeile 46–66
* *GeoMathe_Rechner_Advanced:* Zeile 46–66

---

### 🟨 3. Tab – 1. Hauptaufgabe

* Auswahl des Ausgangspunktes über Dropdown
* Manuelle Eingabe von:

  * Strecke **S**
  * orientierter Richtung **Nu**

**Code‑Referenz:**

* *GeoMathe_Rechner:* Zeile 68–82
* *GeoMathe_Rechner_Advanced:* Zeile 68–82

---

### 🟧 4. Tab – Halbwinkelsatz

Berechnung aller **drei Winkel eines Dreiecks** mithilfe des Halbwinkelsatzes.

**Benötigte Eingaben:**

* Seite **a** → gegenüber von *α*
* Seite **b** → gegenüber von *β*
* Seite **c** → gegenüber von *γ*

**Code‑Referenz:**

* *GeoMathe_Rechner:* Zeile 88–140
* *GeoMathe_Rechner_Advanced:* Zeile 88–140

---

### 🟥 5. Tab – Numerisch stabiler Algorithmus

Berechnung des **Neupunkts Pₙ** (Rückwärtsschnitt) inklusive aller Zwischenergebnisse:

* Hilfsgrößen *a / b*
* *λ (lambda)*, *μ (mu)*
* *s²mn* usw.

**Benötigte Eingaben:**

* Punkte **L**, **M**, **R**
* Winkel **α**, **β**

**Code‑Referenz:**

* *GeoMathe_Rechner:* Zeile 146–210
* *GeoMathe_Rechner_Advanced:* Zeile 146–210

---

## 📘 Hinweise

* Die Berechnungen basieren auf dem **GeoMathe‑1‑Skriptum**
* Der Rest der Bedienung ist weitgehend **selbsterklärend**

---

🎉 **Viel Spaß beim Verwenden des GeoMathe Rechners!**

