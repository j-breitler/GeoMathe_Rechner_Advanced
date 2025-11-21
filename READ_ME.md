------------------------------ GeoMathe-Rechner.READ_ME------------------------------
------------------------------ Author: Julian Breitler ------------------------------
-------------------------------------------------------------------------------------

Der hier beschriebene GeoMathe-Rechner basierd auf Python. 
Falls ihr Python nicht installiert habt bitte installieren.

Es gibt 2 Versionen:

o Version 1 (GeoMathe_Rechner)

	Diese Version funktioniert mit der Python Standard Bibliothek,
	das heißt man muss nur Python installiert haben. Sonst nichts.

o Version 2 (GeoMathe_Rechner_Advanced)

	Diese Version hat ein schöneres User Interface (größere Schrift, cleaner, Abbildungen, usw.), benötigt aber
	eine externe Bibliothek/Package (matplotlib).
	Dieses Package kann auf viele verschiedene Arten installiert werden.
	(wie man python packages installiert bitte Googln oder KI fragen)

Öffnen des Rechners:

	Beide Rechner können auf mehreren Wegen geöffnet werden:
	
	o Direkt das .py file öffnen "Doppelklick"
	o .py File über das Terminal/Command Zeile öffnen
	o .py File in einem Interpreter wie VS Code/ PyCharm/ Spyder/ usw. öffnen und ausführen
	o .exe File ausführen "Doppelklick" (nur für Advanced-Rechner, hierbei wird das benötigte matplotlib Package automatisch 		  installiert, funktioniert nur auf Windows).

	Die einfachste/sicherste Lösung ist "Doppelklick" auf das .py File des einfach GeoMathe_Rechners 
	(sollte bei jedem funktionieren sofern Python installiert ist)  
	
Beide Rechner sind folgender Maßen aufgebaut:

1. Tab 
	
	"Point Management" 

		Punkte eingeben und unter einer "Point ID" speichern. 
		Diese gespeicherten Punkte kann man dann im 2. und 3. Tab (2.HA und 1.HA) wieder aufrufen und einsetzen.
2. Tab
	
	2. Hauptaufgabe berechnen

		Die zuvor eingegebenen Punkte können über ein Drop-Down Menu in die 2.HA 
		eingesetzt werden. Der Rechner führt automatisch die "modulo operation" durch falls eine orientierte Richtung
		>400 gon oder negativ ist. Dabei wir aber immer der original Wert angezeigt und es erscheint ein Hinweis darauf
		das die modulo operation durchgeführt wurde. 
		Wer die Berechnung nachvollziehen will kann im Code nachschauen:
			GeoMathe_Rechner ab Zeile: 46-66 
			GeoMathe_Rechner_Advanced ab Zeile: 46-66 
3. Tab

	1. Hauptaufgabe berechnen

		Ähnlich wie bei der 2.HA kann der Ausgangspunkt über ein Drop-Down Menu ausgewählt
		werden. Die Strecke (S) und die orientierte Richtung (Nu) muss manuell eingetragen werden.
		Wer die Berechnung nachvollziehen will kann im Code nachschauen:
			GeoMathe_Rechner Zeile: 68-82
			GeoMathe_Rechner_Advanced ab Zeile: 68-82
4. Tab
	
	Halbwinkelsatz berechnen 

		Da die Berechnung und das eintippen der Formel für den HWS aufwendig/fehleranfällig und
		nervig ist, rechnet diese Funktion alle 3 Winkel des gegebenen Dreiecks aus. Dafür werden als Input die drei Seiten
		benötigt:
		o a
 		o b
		o c 
		Wobei 
		"alpha" gegenüber von "a" liegt, 
		"beta" gegenüber von "b" und 
		"gamma" gegenüber von "c".

		Wer die Berechnung nachvollziehen will kann im Code nachschauen:
			GeoMathe_Rechner ab Zeile: 88-140
			GeoMathe_Rechner_Advanced ab Zeile: 88-140
5. Tab

	Numerisch-stabilen Algorithmus berechnen
	
		Da die Berechnung und das eintippen der Formel für den Numerisch-stabilen Algo. aufwendig/fehleranfällig und
		nervig ist, rechnet diese Funktion den Neupunkt Pn aus. Dabei werden auch alle Zwischenergebnisse ausgegeben
		(Hilfsgrößen a/b, lambda, mu, s^2mn, usw.)
		Als Input benötigt man die Koordinaten der Punkte
		o L
		o M
		o R 
		(siehe Rückwärtsschnitt im Skriptum), sowie
		o alpha
		o beta 

		Wer die Berechnung nachvollziehen will kann im Code nachschauen:
			GeoMathe_Rechner ab Zeile: 146-210
			GeoMathe_Rechner_Advanced ab Zeile: 146-210

	

Ansonsten sollte alles selbsterklärend sein. Die Berechnungen basieren auf dem Skriptum von GeoMathe-1.



